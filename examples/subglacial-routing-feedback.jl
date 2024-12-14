# This is an example on subglacial water routing.  If this does not mean
# anything to you, you can just view this example with phi being an arbitrary
# DEM.
#
# Additionally to subglacial-routing.jl, this adds a feedback where
# waterflow will generate extra melt due to dissipation of potential energy.

# Load packages
#
# to allow plotting PyPlot must be installed in your global environment
if !@isdefined plotyes
    plotyes = true
end
if plotyes
    @eval using GLMakie
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

open(joinpath(@__DIR__, "data/bed.bin"), "r") do f
     read!(f, bed)
end
open(joinpath(@__DIR__, "data/surface.bin"), "r") do f
     read!(f, surface)
end
open(joinpath(@__DIR__, "data/x.bin"), "r") do f
     read!(f, x)
end
open(joinpath(@__DIR__, "data/y.bin"), "r") do f
     read!(f, y)
end
const ddxx = x[2]-x[1]

# Make hydraulic potential (the "effective DEM")
flotation_fraction = 0.95 # subglacial water pressure as fraction of ice overburden pressure
                          # likely ∈ (0.6. 1.0)
day = 24*60*60 # seconds in a day
rho_w = 1000 # density water
rho_i = 910 # density ice
g = 9.81 # accel. due to gravity
# hydraulic potential (in m H2O)
phi = bed + flotation_fraction * rho_i/rho_w * (surface - bed)


# To track the amount of melted ice, we'll use three cellareas as inputs,
# which in turns leads to three outputs in `area`:
# - water discharge due to source and extra melt
# - water discharge due to extra melt only
# - melt rate in m/s at each cell
cellarea = (fill!(similar(surface), 0.1/day*ddxx*ddxx), # summer conditions: 10cm melt per day over the whole glacier
            fill!(similar(surface), 0.0), # flux due to extra melt; no source here
            fill!(similar(surface), 0.0)) # areal source term due to extra melt; no source here

const phi_ = phi # make const for use in below function

"""
    melting(uparea, dir, ij)

Calculate melt due to dissipation of potential energy given by

-Q*∇phi*g / L,

for phi in [m H2O].  This needs to be scaled from volume per length per time
to volume per area per time.

Returns
- water flux including additional melt [m3/s]
- water flux due to additional melt only [m3/s]
- melt rate at each cell [m/s]

Notes:
- the melt rate is negative in places where water is routed
  out of depressions as the water flows up the hydraulic potential.
  However, freezing will always be less than melt.
- note that the areal source is dependent on the grid size: D8 will
  (almost) always concentrate flow into one cell, thus areal-melt
  will increase.  Note however, that the total volume is correct.
"""
function melting(uparea, ij, dir)
    # constants
    L = 333.55e3 # latent heat of fusion

    dirij = dir[ij]
    Q_total, Q_extra_melt = uparea[1], uparea[2]
    ds = iseven(dirij) ? ddxx : sqrt(2)*ddxx
    phi2 = phi_[ij + WWF.dir2ind(dirij)]
    phi1 = phi_[ij]
    ∇phi = isnan(phi2) ? 0.0 : (phi2 - phi1) / ds
    melt = -Q_total*∇phi * g / L * ds
    return Q_total+melt, Q_extra_melt+melt, melt/ddxx^2
end


# # Route the water
area, slen, dir, nout, nin, sinks, pits, c, bnds  = WWF.waterflows(phi, drain_pits=true,
                                                            cellarea,
                                                            feedback_fn=melting)

# find lakes
phi_filled = WWF.fill_dem(phi, sinks, dir) #, small=1e-6)
lake_depth = phi_filled .- phi
# Plot it
if plotyes
    fig = Figure()
    plt_it!(Axis(fig[1,1]), x, y, phi)
    plt_area!(Axis(fig[1,2]), x, y, area[1], pits)
    plt_area!(Axis(fig[1,1]), x, y, area[2], pits)
    plt_catchments!(Axis(fig[1,1]), x, y, c)
    heatmap!(Axis(fig[1,1]), x, y, lake_depth)
end
