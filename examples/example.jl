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

"An artificial DEM"
function ele(x, y; withpit=false, randfac=0.0)
    out = - (x^2 - 1)^2 - (x^2*y - x - 1)^2 + 6 + 0.1*x + 3*y
    if withpit
        out -= 2*exp(-(x^2 + y^2)*50)
    end
    out += randfac*randn()
    return out<0 ? 0.0 : out
end
dx = 0.01
xs = -1.5:dx:1
ys = -0.5:dx:3.0
dem = ele.(xs, ys', randfac=0.1, withpit=true);
plotyes && WWF.plt.heatmap(xs, ys, dem)

area, slen, dir, nout, nin, pits, c, bnds  = WWF.waterflows(dem, drain_pits=true);

@assert size(dem)==(length(xs), length(ys))
plotyes && WWF.plt.plotit(xs, ys, dem)
plotyes && WWF.plt.plotarea(xs, ys, area, pits)

plotyes && WWF.plt.heatmap(xs, ys, c)

demf = WWF.fill_dem(dem, pits, dir) #, small=1e-6)
plotyes && WWF.plt.heatmap(xs, ys, demf.-dem)
