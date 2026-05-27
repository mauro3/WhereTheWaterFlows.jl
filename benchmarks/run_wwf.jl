using Pkg

Pkg.activate(@__DIR__)
benchdir = @__DIR__
cd(benchdir)
include(joinpath(benchdir, "helpers.jl"))

using WhereTheWaterFlows
using Rasters
using CSV
using DataFrames
import ArchGDAL
using Statistics

#Important: Julia cannot reliably change thread count from inside the script
#set JULIA_NUM_THREADS=4 or JULIA_NUM_THREADS=1 before launching.

function load_real_dem(datadir::String, dataset::String)
    bedfile = ensure_small_dem(datadir)

    bed_raster = replace_missing(Raster(bedfile), NaN32)

    if dataset == "small"
        return Matrix(bed_raster)
    end

    bedfile_large = joinpath(datadir, "swissalti3d_tilelarge.tif")
    if !isfile(bedfile_large)
        println("Creating large DEM via Rasters + ArchGDAL:")
        println(bedfile_large)
        bed_raster_large = resample(bed_raster; size = 4 .* size(bed_raster), method="cubic")
        bed_large = Float32.(Matrix(bed_raster_large))

        geotransform = ArchGDAL.read(bedfile) do ds
            ArchGDAL.getgeotransform(ds)
        end
        proj_wkt = ArchGDAL.read(bedfile) do ds
            ArchGDAL.getproj(ds)
        end

        geotransform_large = copy(geotransform)
        geotransform_large[2] /= 4
        geotransform_large[6] /= 4

        ArchGDAL.create(
            bedfile_large;
            driver = ArchGDAL.getdriver("GTiff"),
            width = size(bed_large, 2),
            height = size(bed_large, 1),
            nbands = 1,
            dtype = Float32,
        ) do ds
            ArchGDAL.setgeotransform!(ds, geotransform_large)
            ArchGDAL.setproj!(ds, proj_wkt)
            band = ArchGDAL.getband(ds, 1)
            ArchGDAL.write!(band, bed_large)
        end

        println("Saved large DEM: ", bedfile_large)
        return Matrix(bed_raster_large)
    end

    return Matrix(replace_missing(Raster(bedfile_large), NaN32))
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
    scriptname = "run_wwf.jl",
    defaults_text = "real/small=6 runs, real/large=3 runs, mock=1 run",
    default_runs = default_runs,
    mock_default_runs = 1,
)

nthreads = Threads.nthreads()
println("Julia threads: ", nthreads)
println("Mode: ", opts.mode)
println("Dataset: ", opts.dataset)
println("Runs: ", opts.runs)


datadir = joinpath(benchdir, "data_raw")
outdir = joinpath(benchdir, "outputs", "wwf")


mkpath(datadir)
mkpath(outdir)

if opts.mode == "mock"
    bed = create_mock_dem(opts.dataset)
else
    bed = load_real_dem(datadir, opts.dataset)
end

npixels = count(.!isnan.(bed))

nruns = opts.runs
runtimes = Float64[]

for i in 1:nruns
    println("WWF run $i / $nruns")

    Base.GC.gc() # run garbage collector so it does not run during bench
    t = @elapsed begin
        global area, slen, dir, nout, nin, sinks, pits, c, bnds, extra
        area, slen, dir, nout, nin, sinks, pits, c, bnds, extra =
            waterflows(bed)
    end

    push!(runtimes, t)
    println("  runtime: ", round(t; digits=3), " s")
end

if opts.mode == "mock"
    warmup_drop = 0
elseif opts.dataset == "large"
    warmup_drop = min(1, length(runtimes) - 1)
else
    warmup_drop = min(2, length(runtimes) - 1)
end
runtimes = runtimes[(warmup_drop + 1):end]

println("Mean runtime: ", round(mean(runtimes); digits=3), " s")

# -------------------------------------------------------
#  DEM stats
# -------------------------------------------------------

(; pits) = waterflows(bed, drain_pits=false)
println("\nDEM stats:")
println("Size: $(size(area))")
println("pits: $(length(pits))")

# -------------------------------------------------------
# Save benchmark CSV
# -------------------------------------------------------

mean_runtime = mean(runtimes)
std_runtime  = std(runtimes)

csvfile = joinpath(outdir, "benchmark_wwf.csv")

newrow = DataFrame(
    method = ["wwf"],
    features = ["waterflows ($(opts.mode), $(opts.dataset))"],
    npixels = [npixels],
    nthreads = [nthreads],
    nruns = [nruns],
    runtime_mean_s = [mean_runtime],
    runtime_std_s = [std_runtime]
)

if isfile(csvfile)

    old = CSV.read(csvfile, DataFrame)

    benchmark = vcat(old, newrow)

else

    benchmark = newrow

end

CSV.write(csvfile, benchmark)

println("Saved benchmark CSV: ", csvfile)



#=
# -------------------------
# Plot DEM
# -------------------------
using CairoMakie

area, slen, dir, nout, nin, sinks, pits, c, bnds, extra =
    waterflows(bed)

x = 1:size(bed, 1)
y = 1:size(bed, 2)
fig0 = Figure()
ax0 = Axis(fig0[1, 1], title = "DEM")
hm0 = heatmap!(ax0, x, y, bed)
Colorbar(fig0[1, 2], hm0, label = "Elevation")
save(joinpath(outdir, "wwf_dem.png"), fig0)

# -------------------------
# Plot accumulation
# -------------------------
fig1 = Figure()
ax1 = Axis(fig1[1, 1], title = "Flow accumulation")
hm1 = plt_area!(ax1, x, y, area)
Colorbar(fig1[1, 2], hm1, label = "log10(area)")
save(joinpath(outdir, "wwf_area.png"), fig1)

# -------------------------
# Plot catchments
# -------------------------
fig2 = Figure()
ax2 = Axis(fig2[1, 1], title = "Catchments")
hm2 = plt_catchments!(ax2, x, y, c)
save(joinpath(outdir, "wwf_catchments.png"), fig2)
=#
