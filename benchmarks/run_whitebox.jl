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
import Whitebox as wbt

function write_dem_geotiff(path::String, dem::Matrix{Float32}; epsg::Int = 2056)
    ArchGDAL.create(
        path;
        driver = ArchGDAL.getdriver("GTiff"),
        width = size(dem, 2),
        height = size(dem, 1),
        nbands = 1,
        dtype = Float32,
    ) do ds
        ArchGDAL.setgeotransform!(ds, [0.0, 1.0, 0.0, 0.0, 0.0, -1.0])
        ArchGDAL.setproj!(ds, ArchGDAL.toWKT(ArchGDAL.importEPSG(epsg)))
        band = ArchGDAL.getband(ds, 1)
        ArchGDAL.write!(band, dem)
    end
end

function create_mock_dem_file(outdir::String, dataset::String)
    dem = create_mock_dem(dataset)

    dem_file = joinpath(outdir, "whitebox_mock_$(dataset)_dem.tif")
    write_dem_geotiff(dem_file, dem)
    return dem_file
end

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
    scriptname = "run_whitebox.jl",
    defaults_text = "real/small=6 runs, real/large=3 runs, mock=1 run",
    default_runs = default_runs,
    mock_default_runs = 1,
)

# -------------------------------------------------------
# Paths
# -------------------------------------------------------

datadir = joinpath(benchdir, "data_raw")
outdir = joinpath(benchdir, "outputs", "whitebox")

mkpath(datadir)
mkpath(outdir)

if opts.mode == "mock"
    demfile = create_mock_dem_file(outdir, opts.dataset)
else
    demfile = ensure_real_dem(datadir, opts.dataset)
end

# -------------------------------------------------------
# Benchmark settings
# -------------------------------------------------------

nthreads = Threads.nthreads()
nruns = opts.runs

println("Julia threads: ", nthreads)
println("Mode: ", opts.mode)
println("Dataset: ", opts.dataset)
println("Runs: ", nruns)

println()
println("Whitebox benchmark on $(opts.dataset) DEM")
println("DEM: ", demfile)

# -------------------------------------------------------
# Run benchmark
# -------------------------------------------------------

csvfile = joinpath(outdir, "benchmark_whitebox.csv")

all_rows = DataFrame()

tmpdir = mktempdir("/tmp"; prefix = "wwf_whitebox_")
atexit(() -> rm(tmpdir; force = true, recursive = true))
conditioned_dem_file = joinpath(tmpdir, "whitebox_$(opts.dataset)_conditioned_dem.tif")
accumulation_file = joinpath(tmpdir, "whitebox_$(opts.dataset)_d8_accumulation.tif")

runtimes = Float64[]

for i in 1:nruns
    println("Whitebox run $i / $nruns")

    rm(conditioned_dem_file; force = true)
    rm(accumulation_file; force = true)

    GC.gc()

    t = @elapsed begin
        wbt.breach_depressions(
            dem = demfile,
            output = conditioned_dem_file,
            fill_pits = false,
        )

        wbt.d8_flow_accumulation(
            i = conditioned_dem_file,
            output = accumulation_file,
            out_type = "cells",
            pntr = false,
        )
    end

    push!(runtimes, t)

    println("  runtime: ", round(t; digits = 3), " s")
    println("\n\n")

end

if opts.mode == "mock"
    warmup_drop = 0
elseif opts.dataset == "large"
    warmup_drop = min(1, length(runtimes) - 1)
else
    warmup_drop = min(2, length(runtimes) - 1)
end
runtimes = runtimes[(warmup_drop + 1):end]

mean_runtime = mean(runtimes)
std_runtime = std(runtimes)

println(
    "Mean runtime: ",
    round(mean_runtime; digits = 3),
    " ± ",
    round(std_runtime; digits = 3),
    " s",
)

dem_raster = Raster(demfile)
npixels = count(.!ismissing.(Matrix(dem_raster)))

row = DataFrame(
    method = ["whitebox"],
    features = ["breach_depressions(fill_pits=false) + d8_flow_accumulation ($(opts.mode), $(opts.dataset))"],
    dataset = [opts.dataset],
    npixels = [npixels],
    nthreads = [nthreads],
    nruns = [nruns],
    runtime_mean_s = [mean_runtime],
    runtime_std_s = [std_runtime],
)

append!(all_rows, row)

#=
# Plot outputs
# -------------------------------------------------------
using CairoMakie

dem = Raster(demfile)
acc = Raster(accumulation_file)

x, y = collect.(dims(dem))

fig = Figure()
ax = Axis(fig[1, 1], title = "Whitebox accumulation ($(opts.dataset))")
hm = heatmap!(ax, x, y, permutedims(log10.(Matrix(acc))))
Colorbar(fig[1, 2], hm, label = "log10(accumulation)")
save(joinpath(outdir, "whitebox_$(opts.dataset)_accumulation.png"), fig)

fig = Figure()
ax = Axis(fig[1, 1], title = "Whitebox conditioned DEM ($(opts.dataset))")
hm = heatmap!(ax, x, y, permutedims(Matrix(Raster(conditioned_dem_file))))
Colorbar(fig[1, 2], hm, label = "Elevation")
save(joinpath(outdir, "whitebox_$(opts.dataset)_conditioned_dem.png"), fig)
=#


# -------------------------------------------------------
# Append to CSV
# -------------------------------------------------------

if isfile(csvfile)
    old = CSV.read(csvfile, DataFrame)
    out = vcat(old, all_rows)
else
    out = all_rows
end

CSV.write(csvfile, out)

println()
println("Saved benchmark CSV: ", csvfile)
