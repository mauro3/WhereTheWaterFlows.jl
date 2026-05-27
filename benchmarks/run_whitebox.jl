using Pkg
Pkg.activate(@__DIR__)
root = normpath(joinpath(@__DIR__, ".."))

using Downloads
using CSV
using DataFrames
using Statistics
using Rasters
using CairoMakie
import ArchGDAL
import Whitebox as wbt

# -------------------------------------------------------
# Paths
# -------------------------------------------------------

datadir = joinpath(root, "data_raw")
outdir = joinpath(root, "outputs", "whitebox")

mkpath(datadir)
mkpath(outdir)

dem_small = joinpath(datadir, "swissalti3d_tile.tif")
dem_large = joinpath(datadir, "swissalti3d_tilelarge.tif")

# -------------------------------------------------------
# Download small DEM if missing
# -------------------------------------------------------

if !isfile(dem_small)
    stac_url = "https://data.geo.admin.ch/api/stac/v1/collections/ch.swisstopo.swissalti3d/items?bbox=7.7,46.5,7.8,46.6&limit=20"

    txt = Downloads.download(stac_url) |> read |> String

    urls = String.(m.match for m in eachmatch(r"https://[^\" ]+\.tif", txt))
    candidates = filter(u -> occursin("_0.5_2056_", u), urls)

    if isempty(candidates)
        error("No swissALTI3D 0.5 m GeoTIFF found.")
    end

    url = first(candidates)

    println("Downloading DEM:")
    println(url)

    Downloads.download(url, dem_small)
end

# -------------------------------------------------------
# Check large DEM
# -------------------------------------------------------

if !isfile(dem_large)
    error("""
    Large DEM not found:

    $dem_large

    Create it first, for example from MATLAB or GDAL.
    """)
end

# -------------------------------------------------------
# Benchmark settings
# -------------------------------------------------------

nruns = 6
nthreads = Threads.nthreads()

println("Julia threads: ", nthreads)

datasets = [
    ("small", dem_small),
    ("large", dem_large),
]

# -------------------------------------------------------
# Run benchmark
# -------------------------------------------------------

csvfile = joinpath(outdir, "benchmark_whitebox.csv")

all_rows = DataFrame()

for (dataset, demfile) in datasets

    println()
    println("Whitebox benchmark on $dataset DEM")
    println("DEM: ", demfile)

    pointer_file = joinpath(outdir, "whitebox_$(dataset)_d8_pointer.tif")
    accumulation_file = joinpath(outdir, "whitebox_$(dataset)_d8_accumulation.tif")

    runtimes = Float64[]

    for i in 1:nruns
        println("Whitebox run $i / $nruns")

        rm(pointer_file; force = true)
        rm(accumulation_file; force = true)

        GC.gc()

        t = @elapsed begin
            wbt.d8_pointer(
                dem = demfile,
                output = pointer_file,
                esri_pntr = false
            )

            wbt.d8_flow_accumulation(
            i = pointer_file,
            output = accumulation_file,
            out_type = "cells",
            pntr = true,
            esri_pntr = false
        )
        end

        push!(runtimes, t)

        println("  runtime: ", round(t; digits = 3), " s")
    end

    mean_runtime = mean(runtimes)
    std_runtime = std(runtimes)

    println(
        "Mean runtime: ",
        round(mean_runtime; digits = 3),
        " ± ",
        round(std_runtime; digits = 3),
        " s"
    )

    dem_raster = Raster(demfile)
    npixels = count(.!ismissing.(Matrix(dem_raster)))

    row = DataFrame(
        method = ["whitebox"],
        features = ["d8_pointer + d8_flow_accumulation"],
        dataset = [dataset],
        npixels = [npixels],
        nthreads = [nthreads],
        nruns = [nruns],
        runtime_mean_s = [mean_runtime],
        runtime_std_s = [std_runtime],
    )

    append!(all_rows, row)

    # Plot outputs
    # -------------------------------------------------------

    dem = Raster(demfile)
    acc = Raster(accumulation_file)
    ptr = Raster(pointer_file)

    x, y = collect.(dims(dem))

    fig = Figure()
    ax = Axis(fig[1, 1], title = "Whitebox accumulation ($dataset)")
    hm = heatmap!(ax, x, y, permutedims(log10.(Matrix(acc))))
    Colorbar(fig[1, 2], hm, label = "log10(accumulation)")
    save(joinpath(outdir, "whitebox_$(dataset)_accumulation.png"), fig)

    fig = Figure()
    ax = Axis(fig[1, 1], title = "Whitebox D8 pointer ($dataset)")
    hm = heatmap!(ax, x, y, permutedims(Matrix(ptr)))
    Colorbar(fig[1, 2], hm, label = "D8 pointer")
    save(joinpath(outdir, "whitebox_$(dataset)_pointer.png"), fig)

end


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