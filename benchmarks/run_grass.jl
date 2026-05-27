using Pkg
Pkg.activate(@__DIR__)
benchdir = @__DIR__
cd(benchdir)
include(joinpath(benchdir, "helpers.jl"))

using CSV
using DataFrames
using Statistics
using Rasters
import ArchGDAL

function ensure_real_dem(datadir::String, dataset::String)
    dem_small = ensure_small_dem(datadir)
    dem_large = joinpath(datadir, "swissalti3d_tilelarge.tif")

    if dataset == "small"
        return dem_small
    end

    if !isfile(dem_large)
        error(
            "Large DEM not found: $dem_large\n" *
            "Create it first by running WWF once in real-large mode, e.g.:\n" *
            "julia --project=benchmarks benchmarks/run_wwf.jl --mode real --dataset large --runs 1"
        )
    end

    return dem_large
end

function derive_default_runs(args::Vector{String})
    if "--runs" in args
        return 6
    end
    for i in 1:(length(args) - 1)
        if args[i] == "--dataset" && lowercase(args[i + 1]) == "large"
            return 3
        end
    end
    return 6
end

default_runs = derive_default_runs(ARGS)

opts = parse_common_args(
    ARGS;
    scriptname = "run_grass.jl",
    defaults_text = "real/small=6 runs, real/large=3 runs",
    default_runs = default_runs,
    mock_default_runs = nothing,
)

opts.mode == "real" || error("run_grass.jl only supports --mode real")

grass_bin = get(ENV, "GRASS_BIN", "grass")
grass_nprocs = tryparse(Int, get(ENV, "GRASS_NPROCS", "1"))
isnothing(grass_nprocs) && error("GRASS_NPROCS must be an integer")
grass_nprocs >= 1 || error("GRASS_NPROCS must be >= 1")

function watershed_supports_nprocs(grass_bin::String, mapset::String)
    cmd = `$grass_bin $mapset --exec r.watershed --help`
    help_txt = read(pipeline(cmd; stderr = stdout), String)
    return occursin("nprocs", help_txt)
end

function run_watershed(grass_bin::String, mapset::String, nprocs::Int; supports_nprocs::Bool)
    if supports_nprocs
        run(`$grass_bin $mapset --exec r.watershed -s elevation=dem accumulation=acc drainage=drain nprocs=$nprocs --overwrite`)
    else
        run(`$grass_bin $mapset --exec r.watershed -s elevation=dem accumulation=acc drainage=drain --overwrite`)
    end
end

datadir = joinpath(benchdir, "data_raw")
outdir = joinpath(benchdir, "outputs", "grass")
csvfile = joinpath(outdir, "benchmark_grass.csv")

mkpath(datadir)
mkpath(outdir)

demfile = ensure_real_dem(datadir, opts.dataset)

println("GRASS benchmark on $(opts.dataset) DEM")
println("DEM: ", demfile)
println("Runs: ", opts.runs)
println("GRASS_BIN: ", grass_bin)
println("GRASS_NPROCS: ", grass_nprocs)

tmpdir = mktempdir("/tmp"; prefix = "wwf_grass_")
atexit(() -> rm(tmpdir; force = true, recursive = true))

grassdb = joinpath(tmpdir, "grassdata")
location = joinpath(grassdb, "benchmark_location")
mapset = joinpath(location, "PERMANENT")

mkpath(grassdb)

println("Preparing temporary GRASS location...")
run(`$grass_bin -c $demfile $location --exec r.in.gdal input=$demfile output=dem --overwrite`)
run(`$grass_bin $mapset --exec g.region raster=dem`)

supports_nprocs = watershed_supports_nprocs(grass_bin, mapset)
if supports_nprocs
    println("r.watershed supports nprocs; using GRASS_NPROCS=", grass_nprocs)
else
    println("r.watershed does not support nprocs on this GRASS version; ignoring GRASS_NPROCS")
end

runtimes = Float64[]

for i in 1:opts.runs
    println("GRASS run $i / $(opts.runs)")
    GC.gc()

    t = @elapsed begin
        run_watershed(grass_bin, mapset, grass_nprocs; supports_nprocs = supports_nprocs)
    end

    push!(runtimes, t)
    println("  runtime: ", round(t; digits = 3), " s")
end

if opts.dataset == "large"
    warmup_drop = min(1, length(runtimes) - 1)
else
    warmup_drop = min(2, length(runtimes) - 1)
end
runtimes = runtimes[(warmup_drop + 1):end]

mean_runtime = mean(runtimes)
std_runtime = std(runtimes)

println("Mean runtime: ", round(mean_runtime; digits = 3), " +- ", round(std_runtime; digits = 3), " s")

dem_raster = Raster(demfile)
npixels = count(.!ismissing.(Matrix(dem_raster)))

newrow = DataFrame(
    method = ["grass"],
    features = [
        supports_nprocs ?
        "r.watershed -s (D8, real, $(opts.dataset), nprocs=$(grass_nprocs))" :
        "r.watershed -s (D8, real, $(opts.dataset), nprocs=unsupported)"
    ],
    dataset = [opts.dataset],
    npixels = [npixels],
    nthreads = [supports_nprocs ? grass_nprocs : missing],
    nruns = [opts.runs],
    runtime_mean_s = [mean_runtime],
    runtime_std_s = [std_runtime],
)

if isfile(csvfile)
    old = CSV.read(csvfile, DataFrame)
    out = vcat(old, newrow)
else
    out = newrow
end

CSV.write(csvfile, out)

println("Saved benchmark CSV: ", csvfile)
