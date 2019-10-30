# WhereTheWaterFlows

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mauro3.github.io/WhereTheWaterFlows.jl/stable) -->
<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mauro3.github.io/WhereTheWaterFlows.jl/dev) -->
[![Build Status](https://travis-ci.org/mauro3/WhereTheWaterFlows.jl.svg?branch=master)](https://travis-ci.org/mauro3/WhereTheWaterFlows.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mauro3/WhereTheWaterFlows.jl?svg=true)](https://ci.appveyor.com/project/mauro3/WhereTheWaterFlows-jl)
[![Codecov](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl)
<!-- [![Coveralls](https://coveralls.io/repos/github/mauro3/WhereTheWaterFlows.jl/badge.svg?branch=master)](https://coveralls.io/github/mauro3/WhereTheWaterFlows.jl?branch=master) -->
<!-- [![Build Status](https://api.cirrus-ci.com/github/mauro3/WhereTheWaterFlows.jl.svg)](https://cirrus-ci.com/github/mauro3/WhereTheWaterFlows.jl) -->

This package implements the D8 flow routing algorithm [1] as well as a
basin-filling algorithm, also by [1]. This allows to calculate water
pathways on a digital elevation model (DEM).

So far little efforts have been made to make this fast or memory
efficient.  The algorithm seems to be of order Q(n) where n is the
number of grid points (provided the number of pits is constant).  For
a 1000x1000 map with 8 pits, it runs in 5s on my laptop from 2012.

## Manual

```julia
using WhereTheWaterFlows, PyPlot
const WWF = WhereTheWaterFlows

"Synthtic DEM with a few maxs and mins"
function peaks2(n=100, randfac=0.05)
    coords = range(-pi, pi, length=n)
    return coords, coords, sin.(coords) .* cos.(coords') .-
        0.7*(sin.(coords.+1) .* cos.(coords')).^8 .+
        randfac*randn(n,n)
end
x,y,dem = peaks2()
area, slen, dir, nout, nin, pits, c, bnds = waterflows(dem)

# log-upslope area as well as pits (sinks)
plotarea(x, y, area, pits)
# catchments
figure()
WWF.heatmap(x,y,c)
# stream length
figure()
WWF.heatmap(x,y,slen)

demf = fill_dem(dem, pits, dir)
# "lake-depth"
figure()
WWF.heatmap(x,y,demf.-dem)
```

# References
[1] O’Callaghan, J. and Mark, D.: The extraction of drainage networks
    from digital elevation data, Comput. Vision Graph., 28, 323–344,
    1984.
