# This is an example of a synthetic valley glacier (a rather big one) with an
# overdeepening.


using WhereTheWaterFlows
const WWFS = WhereTheWaterFlows.Subglacially
const WWF = WhereTheWaterFlows
const WWFR = WhereTheWaterFlows.Randomly
using CairoMakie

function surftopo(x, y; a1=2000/sqrt(1e5), a2=0*8/500^2, s1=200/maximum(x), off=0, xm=x[end÷2],
                  randamp=1)
    xmax, xmin = extrema(x)
    ymax, ymin = extrema(y)
    @assert x[1] == 0
    @assert y[end÷2 + 1] == 0

    rnd = randamp * randn(length(x),length(y))

    # a1 = ice_thick_para.a1;  # factor of sqrt(x) factor
    # s1 = ice_thick_para.s1;  # slope in x-dir of linear part
    # a2 = ice_thick_para.a2;  # factor of parabola a2*y^2
    # xm = ice_thick_para.xm;  # x-coord where it goes from convex to concave

    a1*sqrt.(x) .+ s1*x .+ a2*(x .- xm)/(xmax-xmin).*y.^2 .+ off .+ rnd
end

function bedtopo(x, y; slope=0, y_exp=4, curvature_bed=2000/3.5e3^y_exp,
                 od_xloc1 = 15e3,   # midpoint
                 od_xloc2 = 35e3,   # start point
                 od_yloc1 = -.2e3, # end point
                 od_yloc2 = .2e3,  # start point
                 xsl = 0.2,
                 od_xslope1 = -xsl,
                 od_xslope2 = xsl, # this does not depend on surface slope.
                 od_yslope1 = -0.5,
                 od_yslope2 = 0.5,
                 od_depth = 280, #500;#180;
                 )

    xmax, xmin = extrema(x)
    ymax, ymin = extrema(y)
    @assert x[1] == 0
    @assert y[end÷2 + 1] == 0

    out = slope * x .+ curvature_bed * abs.(y).^y_exp


    # overdeepening:
    if od_depth!=0
        xwidth1 = abs( 0.5*pi/od_xslope1*od_depth);
        startpt = od_xloc1 + xwidth1/2;
        xdir = 1 .- (cosdown.(x, startpt,     od_xslope1/od_depth) + cosup.(x, od_xloc2, od_xslope2/od_depth));
        ydir = 1 .- (cosdown.(y, od_yloc1, od_yslope1/od_depth) + cosup.(y, od_yloc2, od_yslope2/od_depth));

        odeep = xdir.*ydir;

        odeep = - od_depth* odeep;
    else
        odeep = 0;
    end

    return out + odeep
end

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

dx = 50
x = 0:dx:50e3
y = -3.5e3:dx:3.5e3
day = 24*60*60

surfdem = surftopo(x,y') #, randamp=0);
beddem = bedtopo(x,y');
mask = surfdem.>beddem;
source = trues(size(surfdem))*0.01/day;

# Topo
# heatmap(x,y, surfdem)
# heatmap(x,y,beddem)

@time out = WWFS.waterflows_subglacial(surfdem, beddem, dx, gamma=[0,WWFS.GAMMA][2], drain_pits=true, avoid_sc=false)
(;area, sinks, phi) = out.routing
WWF.plt_area(x, y, area.total)
WWF.plt_area(x, y, out.pressmelt.sc_locs)


# @time out = WWFS.waterflows_subglacial(surface, beddem, dx, 1, source, mask, gamma=[0,WWFS.GAMMA][2])
# WWF.plt_area(x, y, area.total; sinks)

@time out = WWF.waterflows(phi, drain_pits=true);
WWF.plt_area(x, y, out.area)

### Stochastic
cov_fn = [WWFR.gaussian_kernel, WWFR.exponential_kernel][1]
surfdem_uc = WWFR.Uncertainty(absuc=0, reluc=0.005, correlation_length=200, covariance_fn=cov_fn)
beddem_uc = WWFR.Uncertainty(absuc=0, reluc=0.05, correlation_length=200, covariance_fn=cov_fn)
floatfrac = 0.97 * ones(size(surfdem));
floatfrac_uc = WWFR.Uncertainty(absuc=0, reluc=0.1, correlation_length=200, covariance_fn=WWFR.gaussian_kernel)
source = ones(size(surfdem));
source_uc = WWFR.Uncertainty()

# two catchment sinks: one all of x==0, one only half of x=0
ctch_sinks = [CartesianIndices((1:10, 1:length(y)))[:], CartesianIndices((1:10, 1:length(y)÷2))[:]]


model, sample, reduce! = WWFR.make_fns_subglacial(step(x),
                                                   surfdem, surfdem_uc,
                                                   beddem, beddem_uc,
                                                   floatfrac, floatfrac_uc,
                                                   source, source_uc,
                                                   ctch_sinks;
                                                   gamma=[0,WWFS.GAMMA][2],
                                                   mask)
input, output = model(sample()...)
aggr = WWFR.map_mc(model, sample, reduce!, 20)
WWF.plt_area(x, y, aggr.areas_total)
# WWF.plt_area(x, y, aggr.sc_locs, prefn=x->x)
WWF.plt_area(x, y, aggr.catchments[:,:,2], prefn=x->x)
