# This example uses a real DEM (around Unteraar glacier in Switzerland)
# and demonstrates plotting

# Load packages
using WhereTheWaterFlows, CairoMakie
const WWF = WhereTheWaterFlows

# Read data for ??? glacier
nx, ny = 539, 459
surface = Matrix{Float32}(undef, nx, ny)
x = Vector{Float32}(undef, nx)
y = Vector{Float32}(undef, ny)

open(joinpath(@__DIR__, "../data/surface.bin"), "r") do f
     read!(f, surface)
end
open(joinpath(@__DIR__, "../data/x.bin"), "r") do f
     read!(f, x)
end
open(joinpath(@__DIR__, "../data/y.bin"), "r") do f
     read!(f, y)
end

# Route the water
(;area, slen, dir, nout, nin, sinks, pits, c, bnds)  = WWF.waterflows(surface, drain_pits=true)

# Plot it
display(plt_area(x, y, area; sinks))

# plot catchments bigger than 20^2 pixels
display(plt_catchments(x, y, c, minsize=20^2))

# plot lakes (or puddles for this DEM)
surface_filled = WWF.fill_dem(surface, sinks, dir)
display(heatmap(x, y, surface_filled .- surface))
