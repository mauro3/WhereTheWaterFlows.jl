# This example shows routing on a synthetic ice sheet margin.
# It showcases:
# - impact of overdeepening on subglacial water flow 
# - stochastic routing with WWFR

using WhereTheWaterFlows
const WWFS = WhereTheWaterFlows.Subglacially
const WWF = WhereTheWaterFlows
const WWFR = WhereTheWaterFlows.Randomly
using CairoMakie

# SHMIP topo
# https://shmip.bitbucket.io/instructions.html
# https://bitbucket.org/shmip/hydro_intercomparison/src/master/input_functions/topography/topos.jl

"""
Ice sheet margin-like topography.

100km x 20km
"""
sqrt_surface(x,y; randamp=0.2) = 6*( sqrt(x+5e3) - sqrt(5e3) ) + 1 + randamp * randn()
sqrt_bed(x,y) = 0
"""
Bed with a trough running at an angle
"""
function sqrt_bed_trough(x,y; depth=100, slope=8e-4, yslope=4e-4, angle=2, ymidpoint_xloc=40e3, bottomwidth=3e3)
    # xlocation of trough middle for the input y
    xloc = (y - 10e3 + angle*ymidpoint_xloc) / angle
    return trough(x - xloc, y, depth, slope, yslope, bottomwidth)
end
trough(x, y, depth, slope, yslope, bottomwidth=0) =
    depth * ( -1
              + cosup(x-bottomwidth/2, 0, slope)
              + cosup(-x-bottomwidth/2, 0, slope)  # into the trough
              ) * cosup(y, 5e3, yslope)

cosdown(x, end_, slope) = coss(-x, -end_, -slope)
cosup(x, start, slope) = coss(x, start, slope)
function coss(x, start, slope)
    # Cos change from 0 to +/-1
    w = 0.5*pi/slope; # width
    end_ = start + abs(w);

    x<start && return 0.0
    x>end_ && return sign(slope)*1.0

    d = 0.5 - 0.5*cos((x-start)*pi/w);
    return sign(slope)*d
end

dx = 100
x = 0:dx:100e3
y = 0:dx:20e3

surfdem = sqrt_surface.(x,y')
beddem = sqrt_bed_trough.(x,y')
nn = 10
ns = length(y)÷nn
ctch_sinks = [CartesianIndex.(2, (i-1)*ns+1:i*ns) for i=1:nn-1] # avoid edge, for some reason
push!(ctch_sinks, CartesianIndex.(2, (nn-1)*ns+1:length(y)) )

heatmap(x,y,beddem)
heatmap(x,y,surfdem)

# Pressure melting point turned on 
@time out = WWFS.waterflows_subglacial(surfdem, beddem, dx, gamma=[0,WWFS.GAMMA][2], bnd_as_sink=true, drain_pits=true)
(;dir, c, pits, area, sinks) = out.routing
(;sc_locs) = out.pressmelt
cc = WWF.catchments(dir, ctch_sinks);
WWF.plt_area(x, y, area.total; sinks)
# WWF.plt_area(x, y, c; sinks, prefn=x->x)
# WWF.plt_area(x, y, sc_locs, prefn=x->x)
WWF.plt_area(x, y, cc; prefn=x->x, title="pressmelt")
cc_size_g = [sum(cc.==n) for n=1:length(ctch_sinks)]

# Pressure melting point turned off 
@time out = WWFS.waterflows_subglacial(surfdem, beddem, dx, gamma=[0,WWFS.GAMMA][1], bnd_as_sink=true, drain_pits=true)
(;dir, pits, area, sinks) = out.routing
(;sc_locs) = out.pressmelt
cc = WWF.catchments(dir, ctch_sinks);
WWF.plt_area(x, y, area.total; sinks)
# WWF.plt_area(x, y, sc_locs, prefn=x->x)
WWF.plt_area(x, y, cc; sinks, prefn=x->x, title="no press-melt")
cc_size_g0 = [sum(cc.==n) for n=1:length(ctch_sinks)]

# plot size of catchments defined by ctch_sinks
p = plot(1:length(ctch_sinks), cc_size_g)
plot!(1:length(ctch_sinks), cc_size_g0) #, title="Catchment size")#, legend=["press-melt", "no press-melt"])
p

### Stochastic
surfdem = sqrt_surface.(x,y', randamp=0)
#beddem = sqrt_bed.(x,y)
beddem = sqrt_bed_trough.(x,y',depth=500)

cov_fn = [WWFR.gaussian_kernel, WWFR.exponential_kernel][2]
surfdem_uc = WWFR.Uncertainty(absuc=0, reluc=0.005, correlation_length=1e3, covariance_fn=cov_fn)
beddem_uc = WWFR.Uncertainty(absuc=0, reluc=0.05, correlation_length=0.5e3, covariance_fn=cov_fn)
floatfrac = 0.95 * ones(size(surfdem));
floatfrac_uc = WWFR.Uncertainty(absuc=0, reluc=0.05, correlation_length=0.5e3, covariance_fn=cov_fn)
basal_melt = ones(size(surfdem));
basal_melt_uc = WWFR.Uncertainty()
mask = trues(size(surfdem))
model, sample, reduce! = WWFR.make_fns_subglacial(step(x),
                                                  surfdem, surfdem_uc,
                                                  beddem, beddem_uc,
                                                  floatfrac, floatfrac_uc,
                                                  basal_melt, basal_melt_uc,
                                                  ctch_sinks;
                                                  mask)
@time aggr = WWFR.map_mc(model, sample, reduce!, 10);
WWF.plt_area(x, y, aggr.areas_total; sinks)
WWF.plt_area(x, y, aggr.lake_depth_free_surface; prefn=x->x)
# WWF.plt_area(x, y, aggr.sc_locs, prefn=x->x)
hist(aggr.catchment_fluxes.total[2])
WWF.plt_area(x, y, aggr.catchments[:,:,2], prefn=x->x)
