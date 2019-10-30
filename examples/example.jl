# Also needs packages Parameters and VAWTools

plotyes = true
using PyPlot
using WhereTheWaterFlows
if !(@isdefined WWF)
    const WWF = WhereTheWaterFlows
end

"A artificial DEM"
function ele(x, y; withpit=false, randfac=0.0)
    out = - (x^2 - 1)^2 - (x^2*y - x - 1)^2 + 6 + 0.1*x + 3*y
    if withpit
        out -= 2*exp(-(x^2 + y^2)*50)
    end
    out += randfac*randn(size(out))
    return out<0 ? 0.0 : out
end
dx = 0.01
xs = -1.5:dx:1
ys = -0.5:dx:3.0
dem = ele.(xs, ys', randfac=0.1, withpit=true);
plotyes && WWF.heatmap(xs,ys,dem)
area, slen, dir, nout, nin, pits  = WWF.waterflows(dem, fillpits=true);
@assert size(dem)==(length(xs), length(ys))
plotyes && WWF.plotit(xs,ys,dem)
plotyes && WWF.plotarea(xs, ys, area, pits)

c, bnds = WWF.catchments(dir, pits);
plotyes && WWF.heatmap(xs, ys, c)

demf = WWF.fill_dem(dem, pits, dir) #, small=1e-6)
plotyes && WWF.heatmap(xs,ys, demf.-dem)
# WWF.plotarea(xs, ys, demf) # will show some artefacts in filled-in area unless `small` is set

## Sorry, this only works on my machine...
import VAWTools
using Parameters
gorner = VAWTools.read_agr("/home/mauro/model_runs/glads/gorner/input/topo_gorner/data/dhm_gorezg2007.grid")
@unpack x,y,v = gorner;
#v = v + rand(size(v)...)*1e-2;
#vv = VAWTools.boxcar(v, 4); #  + rand(size(v)...)*1e-1; # add some noise to kill lakes
dem = v;

area, slen, dir, nout, nin, pits  = WWF.waterflows(dem);
c, bnds = WWF.catchments(dir, pits);
plotyes && WWF.plotarea(x,y,area, pits)

demf = WWF.fill_dem(dem, pits, dir)
plotyes && WWF.heatmap(x,y, demf.-dem)


dir_, nin_, nout_, pits_, c_, bnds_ = WWF.drainpits(dem, dir, nin, nout, pits);
area_, slen_ = WWF._waterflows(dir_, nout_, nin_, pits_);
WWF.plotarea(x, y, area_, pits_)

using BenchmarkTools
@btime WWF.d8dir_feature(v) # 51.695 ms (388933 allocations: 83.48 MiB)
#@btime WWF.d8dir_feature1(v)

if false # small
    @unpack x,y,v = VAWTools.downsample(gorner, 40);
    WWF.plotarea(x,y,v, prefn=x->x);
    dir, nout, nin, pits = WWF.d8dir_feature(v);
    area, slen, dir, nout, nin, pits  = WWF.waterflows(v);
    c, bnds = WWF.catchments(dir, pits);

    figure(); WWF.heatmap(x,y,slen)
    figure(); WWF.heatmap(x,y,nin)
    figure(); WWF.heatmap(x,y,c)
    figure(); WWF.plotdir(x,y,dir)

    dir_, nin_, nout_, pits_, c_, bnds_ = WWF.drainpits(v, dir, nin, nout, pits);
    area_, slen_ = WWF.waterflows(dir_, nout_, nin_, pits_);
    #figure(); plotdir(x,y,dir_); figure(); plotdir(x,y,dir)
    WWF.plotarea(x,y,area_, pits_)
end
