# WhereTheWaterFlows.jl

[![Build Status](https://github.com/mauro3/WhereTheWaterFlows.jl/workflows/CI/badge.svg)](https://github.com/mauro3/WhereTheWaterFlows.jl/actions)
[![Coverage](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl)
[![DOI](https://zenodo.org/badge/218504028.svg)](https://doi.org/10.5281/zenodo.7086860)

Calculate water flow paths on digital elevation models (DEMs) using the D8
algorithm.  The package delineates catchments, computes upslope area, stream
length, and can fill depressions.

## Installation

```julia
using Pkg
Pkg.add("WhereTheWaterFlows")
```

Plotting functions are provided through a Makie extension.  Load any Makie
backend alongside the package to enable them:

```julia
using WhereTheWaterFlows, GLMakie   # or CairoMakie, WGLMakie, …
```

## Quick start

```julia
using WhereTheWaterFlows

# build a small synthetic DEM
n = 100
xs = range(-π, π, length=n)
dem = sin.(xs) .* cos.(xs')

area, slen, dir, nout, nin, sinks, pits, c, bnds = waterflows(dem)
```

See the [Tutorial](@ref) for a full walkthrough.

## References

[1] O'Callaghan, J. and Mark, D.: The extraction of drainage networks from
digital elevation data, Comput. Vision Graph., 28, 323–344, 1984.

[2] Braun, J. and Willett, S. D.: A very efficient O(n), implicit and parallel
method to solve the stream power equation governing fluvial incision and
landscape evolution, Geomorphology, 180–181, 2013.
[doi:10.1016/j.geomorph.2012.10.008](https://doi.org/10.1016/j.geomorph.2012.10.008)
