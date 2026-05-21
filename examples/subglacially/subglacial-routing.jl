# This routes water on Unteraar glacier, Switzerland.

# Load packages
using WhereTheWaterFlows, CairoMakie
const WWF = WhereTheWaterFlows
const WWFS = WhereTheWaterFlows.Subglacially

# Read data for ??? glacier
nx, ny = 539, 459
beddem = Matrix{Float32}(undef, nx, ny)
surfdem = Matrix{Float32}(undef, nx, ny)
x = Vector{Float32}(undef, nx)
y = Vector{Float32}(undef, ny)
open(joinpath(@__DIR__, "../data/bed.bin"), "r") do f
     read!(f, beddem)
end
open(joinpath(@__DIR__, "../data/surface.bin"), "r") do f
     read!(f, surfdem)
end
open(joinpath(@__DIR__, "../data/x.bin"), "r") do f
     read!(f, x)
end
open(joinpath(@__DIR__, "../data/y.bin"), "r") do f
     read!(f, y)
end
dx = x[2]-x[1]

@time out = WWFS.waterflows_subglacial(surfdem, beddem, dx, gamma=[0,WWFS.GAMMA][2], drain_pits=true, avoid_sc=false)
(;area, sinks, phi, c) = out.routing
WWF.plt_area(x, y, area.total)
WWF.plt_area(x, y, out.pressmelt.sc_locs)

# Plot it
plt_catchments(x, y, c)

display(heatmap(x, y, out.lakes.depth_free_surface))
