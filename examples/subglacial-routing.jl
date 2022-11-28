# This is an example on subglacial water routing.  If this does not mean
# anything to you, you can just view this example with phi being an arbitrary
# DEM.

# Load packages
#
# to allow plotting PyPlot must be installed in your global environment
if !@isdefined plotyes
    plotyes = true
end
if plotyes
    @eval using PyPlot
end
using WhereTheWaterFlows
if !(@isdefined WWF)
    const WWF = WhereTheWaterFlows
end


# Read data for ??? glacier
nx, ny = 539, 459
bed = Matrix{Float32}(undef, nx, ny)
surface = Matrix{Float32}(undef, nx, ny)
x = Vector{Float32}(undef, nx)
y = Vector{Float32}(undef, ny)

open("data/bed.bin", "r") do f
     read!(f, bed)
end
open("data/surface.bin", "r") do f
     read!(f, surface)
end
open("data/x.bin", "r") do f
     read!(f, x)
end
open("data/y.bin", "r") do f
     read!(f, y)
end

# Make hydraulic potential (the "effective DEM")
flotation_fraction = 0.95 # subglacial water pressure as fraction of ice overburden pressure
                          # likely ∈ (0.6. 1.0)
rho_w, rho_i = 1000, 910 # density water & ice
g = 9.81 # accel. due to gravity
# hydraulic potential (in meter-H2O)
phi = bed + flotation_fraction * rho_i/rho_w * (surface - bed)

# Route the water
area, slen, dir, nout, nin, pits, c, bnds  = WWF.waterflows(phi, drain_pits=true)

# Plot it
plotyes && WWF.plotit(x, y, phi)
plotyes && WWF.plotarea(x, y, area, pits)

plotyes && WWF.heatmap(x, y, c)

phi_filled = WWF.fill_dem(phi, pits, dir) #, small=1e-6)
plotyes && WWF.heatmap(x, y, phi_filled .- phi)
