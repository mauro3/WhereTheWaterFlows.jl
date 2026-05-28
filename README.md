# WhereTheWaterFlows

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://mauro3.github.io/WhereTheWaterFlows.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://mauro3.github.io/WhereTheWaterFlows.jl/dev)

[![Build Status](https://github.com/mauro3/WhereTheWaterFlows.jl/workflows/CI/badge.svg)](https://github.com/mauro3/WhereTheWaterFlows.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mauro3/WhereTheWaterFlows.jl?svg=true)](https://ci.appveyor.com/project/mauro3/WhereTheWaterFlows-jl)
[![Coverage](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl)
[![DOI](https://zenodo.org/badge/218504028.svg)](https://doi.org/10.5281/zenodo.7086860)

This package calculates water flow paths on digital elevation models (DEMs).

It implements the D8 flow routing algorithm combined with a
breach-type basin-filling algorithm as descibed by [1].  Its implementation uses a
O(n), recursive algorithm [2].

The model's preformace is on par or better than other routing tools according to our
take-them-with-a-grain-of-salt [benchmarks](benchmarks/README.md).


WhereTheWaterFlows has been used in glaciological contexts, see [References](https://mauro3.github.io/WhereTheWaterFlows.jl/dev/references/), but could be useful in other settings as well.

# Example

Example of an upslope area calculation on a synthetic digital elevation model:

```julia
using WhereTheWaterFlows, CairoMakie

# build a small synthetic DEM
n = 200
x = y = range(-π, π, length=n)
dem = sin.(x).*cos.(y') .+ 0.05*rand(n,n)

out = waterflows(dem)
plt_area(x, y, out.area)
```
![Upslope area](https://github.com/user-attachments/assets/5bfd3617-07c0-4d4c-b5d8-ca9e91a4d956)

For details, see the [documentation](https://mauro3.github.io/WhereTheWaterFlows.jl).

# References
[1] O’Callaghan, J. and Mark, D.: The extraction of drainage networks
    from digital elevation data, Comput. Vision Graph., 28, 323–344,
    1984. [download via google scholar](https://scholar.google.ch/scholar?hl=en&as_sdt=0%2C5&q=The+extraction+of+drainage+networks+from+digital+elevation+data&btnG=)

[2] Braun, J. and Willett, S.D.: A very efficient O(n), implicit and
    parallel method to solve the stream power equation governing
    fluvial incision and landscape evolution, Geomorphology, 2013.
    [link](https://doi.org/10.1016/j.geomorph.2012.10.008)
