module Subglacially
using ..WhereTheWaterFlows: WhereTheWaterFlows, PIT, SINK, BARRIER, dir2ind, diagonal_fac, d8dir_feature,
    on_outer_boundary, ind2dir, iterate_D9, flowsinto, diagonal_fac
using StaticArrays
using ProgressMeter
using GeometryBasics: Point2
using LinearAlgebra, Statistics
using StatsBase: mode
const WWF = WhereTheWaterFlows

const RHOI = 910
# note: BedMachine Antarctica uses 917 and uses ice-equivalent surface elevations. Eq S7
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41561-019-0510-8/MediaObjects/41561_2019_510_MOESM1_ESM.pdf

const RHOW = 1000
const LFUSION = 333.55e3 # latent heat of fusion

const GRAV = 9.81
const GAMMA = -0.31 # Röthlisberger constant (corresponding to ct = -7.4 × 10−8 K Pa−1)

"""
    d8dir_pressmelt((phi, phim), bnd_as_sink, nan_as_sink, gamma, avoid_sc)

D8 directions of a DEM and drainage features, taking pressure-melting point effects
into account (i.e. flow deflection and (almost) no flow when the supercooling
threshold is met).

Elevations with NaN map to dir==BARRIER, cells around them will be set to SINK if nan_as_sink==true
or receive no special treatment otherwise.

The argument `bnd_as_sink` determines whether cells at the domain boundary act as sinks.

Places where supercooling occurs can be treated in two ways: water routes around them
(avoid_sc==true, this sets supercooled cells to BARRIER) or water still flows through then (avoid_sc==false).

Args
- gamma -- the Röthlisberger constant.  Set to 0 to get no supercooling, typical value is `GAMMA`.
- avoid_sc -- if set, then they are set as barrier cells

Return
- dir  - direction, encoded as `dirnums`
- nout - number of outflow cells of a cell (0 or 1)
- nin  - number of inflow cells of a cell (0-8)
- sinks - location of sinks as a `Vector{CartesianIndex{2}}` (sorted) (here dir==SINK)
- pits - location of pits as a `Vector{CartesianIndex{2}}` (sorted) (here dir==PIT)
- dem - DEM, unchanged
- flowdir_extra_output -- location of the supercooled cells (more precisely: the outflow edge
                          of such a cell is supercooled)
"""
function d8dir_pressmelt(phi_phim, bnd_as_sink, nan_as_sink, gamma, avoid_sc)
    phi, phim = phi_phim
    # outputs
    dir = fill!(similar(phi, Int8), 0)
    nout = fill!(similar(phi, Bool), false)
    #sc_locs = fill!(similar(phi, Bool), false)
    sc_locs = phi.==nothing # create a suitable BitArray-like by doing a comparison to something abitrary
    nin = fill!(similar(phi, Int8), 0)
    pits = CartesianIndex{2}[]
    kappas = fill!(similar(phi, Int8), 0)
    dir_ = fill!(similar(phi, Int8), 0) # dir before deflection

    R = CartesianIndices(size(phi))
    Iend = last(R)

    # if there is no supercooling, be done with it
    if gamma==0
        dir, nout, nin, sinks, pits, dem = d8dir_feature(phi, bnd_as_sink, nan_as_sink)
        return dir, nout, nin, sinks, pits, dem, (sc_locs, kappas, dir)
    end

    sc_threshold = (1+gamma) / gamma

    ## Figure out steepest descent of both phi and phim including deflection
    ## and (potentially) no-flow on supercooled cells

    ## Step 1: Get dir for phi and phim
    dir_phi  = d8dir_feature(phi,  bnd_as_sink, nan_as_sink)[1]
    dir_phim = d8dir_feature(phim, bnd_as_sink, nan_as_sink)[1]

    ## Step 2: mark cells with supercooled outflow edge as BARRIER if avoid_sc==true.
    # this can only happen when phi and phim are co-linear and opposite,
    # otherwise flow will just be deflected.
    extra_barriers = CartesianIndex{2}[]
    for I in R
        d, dm = dir_phi[I], dir_phim[I]
        (d==PIT || dm==PIT) && continue
        d>=SINK && continue
        i, im = dir2ind(d), dir2ind(dm)
        if i+im==CartesianIndex(0,0) # co-linear and opposite
            ratio = abs(phim[I+im] - phim[I])*diagonal_fac(im) / (abs(phi[I+i] - phi[I]) *  diagonal_fac(i))
            if ratio > sc_threshold
                sc_locs[I] = true
                if avoid_sc
                    push!(extra_barriers, I)
                end
            end
        end
    end
    if length(extra_barriers)>0
        # update dir_phi
        dir_phi = d8dir_feature(phi, bnd_as_sink, nan_as_sink, CartesianIndex{2}[], extra_barriers)[1]
    end

    # make dir_phi default direction
    dir = copy(dir_phi)

    # Step 3:
    # get dir for all points:
    # - mark supercooled cells (maybe make the barriers)
    # - make the press-melt dir from the two
    # - set dir to dir_phi or to deflected dir
    for I in R
        phidir, phimdir = dir_phi[I], dir_phim[I]
        # no bed-slope means no deflection
        phimdir==PIT && continue
        phimdir==SINK && continue
        phidir==BARRIER && continue

        # calculate angles and gradient ratio
        i, im = dir2ind(phidir), dir2ind(phimdir)
        c, cm = i.I, im.I
        # lambda: angle between -∇ϕ and Q (i.e. `dir`)
        # kappa: angle between -∇ϕ and -∇ϕₘ
        # in units of π/4, instead of dot-product formula use atan2 to also get sign
        lambda = round(Int, (atan(cm[2],cm[1]) - atan(c[2],c[1])) / (pi/4))
        lambda = mod(lambda, -3:4)
        # ratio between gradients, i.e. ratio of difference in respective steepest descent:
        ratio = abs(phim[I+im] - phim[I])*diagonal_fac(im) / (abs(phi[I+i] - phi[I]) * diagonal_fac(i))
        # calculate deflection angle
        kappa = getkappa_D8(ratio, lambda, gamma)
        kappas[I] = kappa

        # adjust the flow direction, i.e. rotate the phidir by kappa
        i_new = i # for abs(kappa)==0
        if 0<abs(kappa)<3
            if abs(kappa)==2
                sinkappa = 1.0 * sign(kappa)
                coskappa = 0.0
                # rotate
                i_new = CartesianIndex(round.(Int, SMatrix{2,2,Float64}(coskappa,  sinkappa,
                                                                        -sinkappa, coskappa) * SVector{2,Float64}(c[1], c[2]))...)
            end
            # if there is flow uphill or into BARRIER with kappa==2, then try kappa==1
            if abs(kappa)==1 || phi[I+i_new]>=phi[I] || dir[I+i_new]==BARRIER
                sinkappa = sqrt(2)/2 * sign(kappa)
                coskappa = sqrt(2)/2
                # rotate

                i_new = CartesianIndex(round.(Int, SMatrix{2,2,Float64}(coskappa,  sinkappa,
                                                                        -sinkappa, coskappa) * SVector{2,Float64}(c[1], c[2]))...)
            end
            # flow uphill or into barrier, thus put it back to no deflection
            if phi[I+i_new]>=phi[I] || dir[I+i_new]==BARRIER
                i_new = i
            end
            if dir[I+i_new]>SINK # revert flow into a BARRIER cell
                i_new = i
            end
        elseif kappa!=0
            error("Something is wrong with the algorithm as abs(kappa)<=2, but got kappa=$kappa.")
        end

        dir_[I] = ind2dir(i)
        dir[I] = ind2dir(i_new)
    end

    return dir, WWF.make_flowfeatures(dir)..., phi, (sc_locs, kappas, dir_)
end

###################################################################
## TO MAKE LOOKUP TABLE OF getkappa_D8 the following code was used:
# import Optim#, Roots
# """
# Function to be maximised to get angle with max melt.
# """
# function L(kappa, R, lambda, gamma)
#     beta = [3/2, 2][1] # turbulent==3/2, laminar==2
#     if kappa<0 || kappa>pi/2
#         return NaN
#     end
#     cos(kappa)^(beta-1) * ( (1+gamma)*cos(kappa) - gamma*cos(lambda-kappa)*R )
# end

# """
# Return kappa (angle between -∇ϕ and Q) which maximizes melt.

# """
# function getkappa(R, lambda, gamma)
#     # do the calculation
#     fac = if lambda<0
#         lambda = -lambda
#         -1
#     else
#         1
#     end

#     sol = Optim.optimize(kappa -> -L(kappa[1], R, lambda, gamma), [1e-3], Optim.LBFGS())
#     if Optim.converged(sol)
#         return fac * Optim.minimizer(sol)[1]
#     else
#         return NaN
#     end
# end
# ##to look at lookup table
# ratio =0:0.01:10
# for l=1:3
#     kappa_ = WWF.getkappa.(ratio, l*pi/4, GAMMA)
#     kappa = round.(Int, kappa_ / (pi/4))
#     plot(ratio, kappa)
#     plot(ratio, kappa_/(pi/4))
#     xlabel("R")
#     ylabel("κ")
# end
# legend(["λ = π/4", "λ = π/2", "λ = 3π/4"])
#
# """Ratio of melt in a "wonkey" channel vs one on a flat bed (i.e. co-linear with -∇ϕ)"""
# meltratio(kappa_max, R, lambda, gamma) = L(kappa_max, R, lambda) / (1+gamma)
# meltratio(R, lambda, gamma) = meltratio(kappa(R, lambda), R, lambda, gamma)
#
# """
# Finds locus where the meltratio==mu as a function of R.
# """
# function locus(mu, R)
#     bracket = (0.0,pi)
#     f = lambda -> meltratio(R, lambda) - mu
#     if sign(f(bracket[1])) == sign(f(bracket[2]))
#         return NaN
#     else
#         return Roots.find_zero(f, bracket)
#     end
# end
#
###############################################################

"""
    getkappa_D8(R, lambda::Int, gamma)

Return kappa (angle between -∇ϕ and Q) which maximizes melt for a given ratio
R = ϕₘ/ϕ. A crude approximation which is exact for D8.

Notes:
- kappa ∈ {-π/2, -π/4, 0, π/4, π/2}
- Only implemented for gamma ∈ {0, -0.31}.
- Takes lambda, the angle between surface and bed slope, in units of π/4, and
  returns kappa in those units too.
- It important is to check that water is not flowing phi-uphill after the deflection.
  If it is the case the deflection needs to be reduced/removed in the calling
  function.
"""
function getkappa_D8(R, lambda::Int, gamma)
    sig = sign(lambda)
    lambda = abs(lambda)
    if gamma==0
        return 0
    elseif gamma==-0.31
        if lambda==0 || lambda==4
            return 0
        elseif lambda==1
            if R<6.677
                return 0
            else
                return sig*1
            end
        elseif lambda==2
            if R<1.512
                return 0
            else
                return sig*1
            end
        elseif lambda==3
            if R<1.274
                return 0
            elseif R<6.663
                return sig*1
            else
                return sig*2
            end
        else
            error("Need 0<=abs(lambda)<=4 but got lambda=$lambda")
        end
    else
        error("not implemented for gamma=$gamma")
    end
end

"""
    phi_fn(surface, thick, floatfrac=1; rhow=RHOW, rhoi=RHOI)

Shreve hydraulic potential [m H2O].
"""
phi_fn(surface, thick, floatfrac=1; rhow=RHOW, rhoi=RHOI) =
    @. floatfrac*thick*rhoi/rhow + (surface - thick)
# Note: there was a waterdepth kwarg
#    @. floatfrac*thick*rhoi/rhow + (surface - thick + waterdepth)
# Maybe there is some value for this?

"""
    lake_depth(phi_filled, phi_orig; fixed_surface=true, rhow=RHOW, rhoi=RHOI)

The lake depth can be calculated assuming
1) adjusting the ice surface (fixed_surface=false)
2) adjusting the ice thickness only but keeping the ice surface fixed (fixed_surface=true)

The formulas are, respectively:

(phi_filled - phi_original)

(phi_filled - phi_original) * (rhow / (rho_w - rho_i) )
"""
function lake_depth(phi_filled, phi_orig; fixed_surface=true, rhow=RHOW, rhoi=RHOI)
    if fixed_surface
        return (phi_filled.-phi_orig).* (rhow / (rhow-rhoi))
    else
        return (phi_filled.-phi_orig)
    end
end

"""
    waterflows_subglacial(surfdem, beddem, dx, 
                               floatfrac=1,
                               source=ones(size(surfdem)),
                               mask=trues(size(surfdem));
                               gamma=-0.31,
                               ctch_sinks=[],
                               rhow = 1000.0,
                               rhoi = 910.0,
                               drain_pits=true,
                               bnd_as_sink=true,
                               nan_as_sink=true)

Does the water flow routing according the D8 algorithm for a subglacial setting using
the Shreve-potential for routing.  Utilizes `WhereTheWaterFlows.waterflows`.

args:
- surfdem, beddem -- the surface and bed DEM
- floatfrac -- flotation fraction
- dx -- grid size (must be equal in x and y-direction)
- source=ones(size(dem)) -- the source per cell, defaults to 1.  If using physical units
                            then use a source in volume per unit area, e.g. m/s
- mask -- routing mask:  where to do routing

kwargs:
- drain_pits -- whether to route through pits (true)
- bnd_as_sink (true) -- whether the domain boundary be sinks, i.e.
                i.e. adjacent cells can drain into them
                or whether to ignore them.
- nan_as_sink (true) -- whether NaNs are sinks
- gamma -- Röthlisberger constant (-0.31, note it's negative)
- avoid_sc (false) -- If set to true, then supercooled cells become BARRIERs.  Note that this then also
                      means that source in those cells will just vanish, i.e. mass is lost.  Probably
                      the default ==false is more physical as waterflow should continue just
                      not in R-channel. Thus recommended to set to false.
- rhow, rhoi -- mean density of water and glacier ice (1000, 910)
- ctch_sinks -- List of sink-areas for which catchments will be calculated.  Each sink can be given as
                a vector of `CartesianIndex` or as `CartesianIndices`.  Note, sink-areas and catchments can be
                overlapping, if desired.

Returns one nested NamedTuple with keys:
- `routing`: fields `area, slen, dir, nout, nin, sinks, pits, c, bnds, phi`
  - `routing.area`: named tuple with `total`, `extra`, `dissipation_melt_rate`, `pressure_melt_rate`
- `pressmelt`: fields `sc_locs, kappas, dir_og`
- `lakes`: fields `depth_fixed_surface, depth_free_surface`
- `sink_catchments`: fields `masks, fluxes`
  - `sink_catchments.fluxes`: named tuple with `total`, `dissipation`, `pressmelt`
"""
function waterflows_subglacial(surfdem::AbstractMatrix, beddem::AbstractMatrix, dx::Number,
                               floatfrac=1,
                               source::AbstractMatrix=fill!(similar(surfdem), 1),
                               mask::AbstractMatrix=fill!(similar(surfdem, Bool), true); # that does not make a BitArray
                               ctch_sinks=[],
                               gamma=GAMMA,
                               rhow=RHOW,
                               rhoi=RHOI,
                               avoid_sc=false,
                               drain_pits=true,
                               bnd_as_sink=true,
                               nan_as_sink=true)
    @assert size(surfdem)==size(beddem)==size(source) "Input arrays not of same size"
    @assert size(floatfrac)==() || size(floatfrac)==size(surfdem)  "Input arrays not of same size"

    # preparation
    thick_orig = (surfdem .- beddem)
    phi  = phi_fn(surfdem, thick_orig, floatfrac; rhow, rhoi)
    phim = phi_fn(surfdem, thick_orig, 0.0,     ; rhow, rhoi)
    phi[.!mask]  .= NaN
    phim[.!mask] .= NaN

    # prepare `cellarea` such that there are four output fields:
    # total flux, extra flux, dissipation-melt-rate, pressure-melt-rate
    inputs = (source*dx^2, fill!(similar(surfdem), 0), fill!(similar(surfdem), 0), fill!(similar(surfdem), 0))

    feedback_fn = make_feedback_fn(phi, phim, gamma, dx)

    area, slen, dir, nout, nin, sinks, pits, c, bnds, (sc_locs, kappas, dir_) =
        WWF.waterflows((phi, phim), inputs;
                       flowdir_fn=(phi_phim, bnd_as_sink, nan_as_sink, _, _) -> d8dir_pressmelt(phi_phim, bnd_as_sink, nan_as_sink, gamma, avoid_sc),
                       drain_pits, bnd_as_sink, nan_as_sink, feedback_fn);


    # dem filling and lakes
    phi_filled = WWF.fill_dem(phi, sinks, dir)
    phi_filled[.!mask] .= NaN;
    lakes = lake_depth(phi_filled, phi; rhow, rhoi)
    lakes_free_surf = lake_depth(phi_filled, phi; fixed_surface=false, rhow, rhoi)

    # Sink areas and fluxes
    #any(isempty.(ctch_sinks)) && error("One of the sinks contains no cells!")
    cs = [WWF.catchment(dir, s) for s in ctch_sinks]
    fluxes = (total = [WWF.catchment_flux((source .+ area[3] .+ area[4])*dx^2, c) for c in cs],
              dissipation = [WWF.catchment_flux(area[3]*dx^2, c) for c in cs],
              pressmelt = [WWF.catchment_flux(area[4]*dx^2, c) for c in cs])

    return (
        routing = (
            area = (
                total = area[1],
                extra = area[2],
                dissipation_melt_rate = area[3],
                pressure_melt_rate = area[4],
            ),
            slen,
            dir,
            nout,
            nin,
            sinks,
            pits,
            c,
            bnds,
            phi,
        ),
        pressmelt = (
            sc_locs,
            kappas,
            dir_og = dir_,
        ),
        lakes = (
            depth_fixed_surface = lakes,
            depth_free_surface = lakes_free_surf,
        ),
        sink_catchments = (
            masks = cs,
            fluxes,
        ),
    )
end

"""
   melt_rates(Q::Number, phi::AbstractMatrix, phim::AbstractMatrix, gamma, ij::CartesianIndex, dir::AbstractMatrix, dx::Number, rhow::Number)

Calculate melt rates due to dissipation of potential energy and due to
pressure melting point effects given by

    (-Q ∇ϕ - γ Q ∇pw)/(GRAV * L),

for ϕ and pw given in [m H2O].

This is calculated for the water flow from cell ij to its downstream cell.

This is used as part of the `feedback_fn` in WWF.waterflows.

Returns
- dissipation melt rate (first term) [m/s]
- pressure melt rate (first term) [m/s]

Note: the dissipation melt rate is negative in places where water is routed
out of depressions as the water flows up the hydraulic potential (this is due
to the breach type algorithm in `drainpits!`).
However, freezing will balance the extra melt occurring when the water descends into
the filled local minimum (pit).
"""
function melt_rates(Q::Number, phi::AbstractMatrix, phim::AbstractMatrix, gamma, ij::CartesianIndex, dir::AbstractMatrix, dx::Number)
    dirij = dir[ij]
    ds = dx/diagonal_fac(dirij)

    phi2 = phi[ij + WWF.dir2ind(dirij)]
    if isnan(phi2)
        ∇phi = 0.0
        ∇pw = 0.0
    else
        phi1 = phi[ij]
        ∇phi = (phi2 - phi1) / ds

        phim2 = phim[ij + WWF.dir2ind(dirij)]
        phim1 = phim[ij]
        ∇pw = ∇phi - (phim2 - phim1) / ds
    end

    disspation_melt_rate = -Q*∇phi * GRAV / LFUSION
    press_melt_rate = -gamma*Q*∇pw * GRAV / LFUSION
    # as Q is in [m3/s], the calculated melt rates are in [m2/s] but we need
    # want to return [m/s] (matching `source` input units in waterflows_subglacial)
    return disspation_melt_rate * ds/dx^2, press_melt_rate * ds/dx^2
end

"""
    make_feedback_fn(phi, phim, gamma, dx)

Make the feedback_fn as used by WWF with the appropriate closures and
taking care of https://github.com/JuliaLang/julia/issues/15276 (if needed)

It returns:
- total discharge
- extra discharge due to dissipation and pressure melting
- dissipation melt rate (m/s)
- pressure melt rate (m/s)
"""
function make_feedback_fn(phi, phim, gamma, dx)
    return function (uparea, ij, dir)
        Q_total, Q_extra = uparea[1], uparea[2]
        Q_total = max(Q_total, 0)
        # don't do  Q_extra = max(Q_extra, 0) as the extra may well be negative
        disspation_melt_rate, press_melt_rate =
            melt_rates(Q_total, phi, phim, gamma, ij, dir, dx)
        return Q_total + disspation_melt_rate*dx^2 + press_melt_rate*dx^2, Q_extra + disspation_melt_rate*dx^2 + press_melt_rate*dx^2, disspation_melt_rate, press_melt_rate
    end
end

include("helpers.jl")
end
