using WhereTheWaterFlows, CairoMakie
const WWF = WhereTheWaterFlows

"An artificial DEM"
function ele(x, y; withlake=false, randfac=0.0)
    out = - (x^2 - 1)^2 - (x^2*y - x - 1)^2 + 6 + 0.1*x + 3*y
    if withlake
        out -= 2*exp(-(x^2 + y^2)*50)
    end
    out += randfac*randn()
    return out<0 ? 0.0 : out
end
dx = 0.01
x = -1.5:dx:1
y = -0.5:dx:3.0
dem = ele.(x, y', randfac=0.1, withlake=true);
display(heatmap(x, y, dem))

(;area, slen, dir, nout, nin, sinks, pits, c, bnds)  = WWF.waterflows(dem, drain_pits=true);

@assert size(dem)==(length(x), length(y))
display(plt_area(x, y, area; sinks))

display(plt_catchments(x, y, c))

demf = WWF.fill_dem(dem, sinks, dir) #, small=1e-6)
display(heatmap(x, y, demf.-dem))
