# WhereTheWaterFlows

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mauro3.github.io/WhereTheWaterFlows.jl/stable) -->
<!-- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mauro3.github.io/WhereTheWaterFlows.jl/dev) -->
[![Build Status](https://travis-ci.com/mauro3/WhereTheWaterFlows.jl.svg?branch=master)](https://travis-ci.com/mauro3/WhereTheWaterFlows.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mauro3/WhereTheWaterFlows.jl?svg=true)](https://ci.appveyor.com/project/mauro3/WhereTheWaterFlows-jl)
<!-- [![Build Status](https://api.cirrus-ci.com/github/mauro3/WhereTheWaterFlows.jl.svg)](https://cirrus-ci.com/github/mauro3/WhereTheWaterFlows.jl) -->
[![Codecov](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl)
<!-- [![Coveralls](https://coveralls.io/repos/github/mauro3/WhereTheWaterFlows.jl/badge.svg?branch=master)](https://coveralls.io/github/mauro3/WhereTheWaterFlows.jl?branch=master) -->


This package implements the D8 flow routing algorithm [1] as well as a
basin-filling algorithm, also by [1].  It uses a O(n), recursive algorithm
similar to [2]. This allows to calculate water
pathways on a digital elevation model (DEM).

This code is reasonably fast: flow routing on a DEM of Antarctica of
about 2e8 points and with 150000 depressions takes about 30s on my
laptop (Ryzen 4750U).

![Upslope area](https://user-images.githubusercontent.com/4098145/67853636-e319b880-fb06-11e9-933d-9f55ace99ce1.png)

Example of upslope area calculated in below example.

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
x,y,dem = peaks2(200)
area, slen, dir, nout, nin, pits, c, bnds = waterflows(dem)

# log-upslope area as well as pits (sinks)
WWF.plotarea(x, y, area, pits)

# log-upslope area over contours of the dem
WWF.plotarea_dem(x, y, dem, area, pits)

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
    1984. [link](http://www.academia.edu/download/47710594/s0734-189x_2884_2980011-020160801-5103-19m2b12.pdf)
[2] Braun, J. and Willett, S.D.: A very efficient O(n), implicit and
    parallel method to solve the stream power equation governing
    fluvial incision and landscape evolution Author links open overlay
    panel [doi](https://doi.org/10.1016/j.geomorph.2012.10.008)
