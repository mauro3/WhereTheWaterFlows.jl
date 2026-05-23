# WhereTheWaterFlows.jl

[![Build Status](https://github.com/mauro3/WhereTheWaterFlows.jl/workflows/CI/badge.svg)](https://github.com/mauro3/WhereTheWaterFlows.jl/actions)
[![Coverage](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl)
[![DOI](https://zenodo.org/badge/218504028.svg)](https://doi.org/10.5281/zenodo.7086860)

WhereTheWaterFlows routes water on gridded topography (or hydraulic potential)
with a D8-style flow model. It supports deterministic flow routing, coupled
feedbacks, subglacial routing physics, and uncertainty propagation.

The package currently has three modules:

- `WhereTheWaterFlows`: core deterministic routing and post-processing (WWF)
- `WhereTheWaterFlows.Subglacially`: subglacial hydraulic-potential routing (WWFS)
- `WhereTheWaterFlows.Randomly`: Monte Carlo wrappers for uncertainty studies (WWFR)

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

Submodules are bundled in the same package and could be used like so:

```julia
using WhereTheWaterFlows
const WWF  = WhereTheWaterFlows
const WWFS = WhereTheWaterFlows.Subglacially
const WWFR = WhereTheWaterFlows.Randomly
```

## Quick start

```julia
using WhereTheWaterFlows

# build a small synthetic DEM
xs = range(-π, π, length=200)
dem = sin.(xs) .* cos.(xs')

out = waterflows(dem)
maximum(out.area), length(out.sinks)
```

See the [Tutorial](@ref) for a full walkthrough of the core API.

## Algorithm

The core routing uses the D8 algorithm: each cell drains to whichever of its
eight neighbours has the steepest downward gradient.  Local minima (pits) are
handled by default via a breach-type algorithm that finds the lowest spillway
for each pit and reverses flow along that path, so the input DEM does not need
to be pre-filled. See O’Callaghan & Mark (1984).

The tree traversal used for accumulation is O(n) and recursive (Braun & Willett, 2013).
On very large DEMs the recursion depth can exceed the default Julia call-stack size and cause
a `StackOverflowError`.  See the `stacksize` keyword argument of `waterflows`
if this occurs.

## Performance

Routing a 14 000 × 14 000 Antarctica DEM (~2 × 10⁸ cells, ~150 000
depressions) takes roughly 30 s on a laptop-class CPU.

## Guides

- [Feedback Functionality](@ref FeedbackGuide)
- [Randomly](@ref RandomlyGuide)
- [Subglacially](@ref SubglaciallyGuide)
- [Examples](@ref)
- [API Reference](@ref)
