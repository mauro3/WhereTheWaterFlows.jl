# WhereTheWaterFlows.jl

[![Build Status](https://github.com/mauro3/WhereTheWaterFlows.jl/workflows/CI/badge.svg)](https://github.com/mauro3/WhereTheWaterFlows.jl/actions)
[![Coverage](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl)
[![DOI](https://zenodo.org/badge/218504028.svg)](https://doi.org/10.5281/zenodo.7086860)

WhereTheWaterFlows routes water on gridded topography (or hydraulic potential)
with a D8-style flow model. It supports deterministic flow routing, coupled
feedbacks, subglacial routing physics, and uncertainty propagation.

The package currently has three modules:

- `WhereTheWaterFlows`: core deterministic routing and post-processing
- `WhereTheWaterFlows.Subglacially`: subglacial hydraulic-potential routing
- `WhereTheWaterFlows.Randomly`: Monte Carlo wrappers for uncertainty studies

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

Submodules are bundled in the same package:

```julia
using WhereTheWaterFlows
WWFS = WhereTheWaterFlows.Subglacially
WWFR = WhereTheWaterFlows.Randomly
```

## Quick start

```julia
using WhereTheWaterFlows

# build a small synthetic DEM
n = 200
xs = range(-π, π, length=n)
dem = sin.(xs) .* cos.(xs')

out = waterflows(dem)
maximum(out.area), length(out.sinks)
```

See the [Tutorial](@ref) for a full walkthrough of the core API.

## Guides

- [Feedback Functionality](@ref FeedbackGuide)
- [Randomly](@ref RandomlyGuide)
- [Subglacially](@ref SubglaciallyGuide)
- [Examples](@ref)
- [API Reference](@ref)
