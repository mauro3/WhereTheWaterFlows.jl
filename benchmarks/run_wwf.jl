using Pkg

Pkg.activate(@__DIR__)
root = normpath(joinpath(@__DIR__, ".."))

using Downloads
using WhereTheWaterFlows
using Rasters
using CairoMakie
using CSV
using DataFrames
import ArchGDAL
using Statistics

#Important: Julia cannot reliably change thread count from inside the script
#set JULIA_NUM_THREADS=4 or JULIA_NUM_THREADS=1 before launching.

nthreads = Threads.nthreads()
println("Julia threads: ", nthreads)


datadir = joinpath(root, "data_raw")
outdir = joinpath(root, "outputs", "wwf")


mkpath(datadir)
mkpath(outdir)

bedfile = joinpath(datadir, "swissalti3d_tile.tif")

if !isfile(bedfile)
    stac_url = "https://data.geo.admin.ch/api/stac/v1/collections/ch.swisstopo.swissalti3d/items?bbox=7.7,46.5,7.8,46.6&limit=20"

    txt = Downloads.download(stac_url) |> read |> String

    urls = String.(m.match for m in eachmatch(r"https://[^\" ]+\.tif", txt))
    candidates = filter(u -> occursin("_0.5_2056_", u), urls)

    url = first(candidates)

    println("Downloading DEM:")
    println(url)

    Downloads.download(url, bedfile)
end

bed_raster = replace_missing(Raster(bedfile), NaN32)
bedfile_large = joinpath(datadir, "swissalti3d_tilelarge.tif")
if ~isfile(bedfile_large)
    bed_raster_large = resample(bed_raster; size=4 .*size(bed_raster))
    # Saving GeoTiff seems to be broken on Julia 1.12 because of
    # https://github.com/JuliaImages/OpenCV.jl/issues/61
    # save(bedfile_large, bed_raster_large)
    #
    # The GDAL command instead is:
    # gdalwarp -ts 8000 8000 -r cubic -co COMPRESS=DEFLATE  swissalti3d_tile.tif swissalti3d_tilelarge.tif
else
    bed_raster_large = replace_missing(Raster(bedfile_large), NaN32)
end
# note: in matlab it is also Float32
bed_small = Matrix(bed_raster)
bed_large = Matrix(bed_raster_large)
bed = [bed_small, bed_large][1]

npixels = count(.!isnan.(bed))

nruns = 6
runtimes = Float64[]

for i in 1:nruns
    println("WWF run $i / $nruns")

    Base.GC.gc() # run garbage collector so it does not run during bench
    t = @elapsed begin
        global area, slen, dir, nout, nin, sinks, pits, c, bnds, extra =
            waterflows(bed)
    end

    push!(runtimes, t)
    println("  runtime: ", round(t; digits=3), " s")
end

runtimes = runtimes[3:end]

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
    features = ["waterflows"],
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



# -------------------------
# Plot DEM
# -------------------------

area, slen, dir, nout, nin, sinks, pits, c, bnds, extra =
    waterflows(bed)

x, y = [d.val.data .- d[1] for d in dims(bed_raster)]
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

