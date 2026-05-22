# [Subglacially](@id SubglaciallyGuide)

`WhereTheWaterFlows.Subglacially` extends core routing for subglacial settings.

Compared with plain `waterflows`, it provides:

- routing based on hydraulic potential,
- pressure-melting/supercooling-aware deflection,
- lake-depth diagnostics,
- optional sink-catchment flux aggregation.

This module is aimed at routing water under ice using Shreve-style hydraulic
potential and optional pressure-melting-point effects.

## Minimal run

```@example subglacially
using WhereTheWaterFlows
using Random

Random.seed!(42)

WWFS = WhereTheWaterFlows.Subglacially

n = 90
dx = 100.0
x = range(0, step=dx, length=n)
X = repeat(collect(x), 1, n)
Y = repeat(collect(x)', n, 1)

surfdem = 1200.0 .+ 0.015 .* X .+ 0.015 .* Y .+ 4.0 .* randn(n, n)
beddem = 950.0 .+ 0.01 .* X .+ 0.01 .* Y .- 100.0 .* exp.(-((X .- x[end] / 2).^2 .+ (Y .- x[end] / 2).^2) ./ (2 * 4000.0^2))
surfdem = max.(surfdem, beddem .+ 10.0)

out = WWFS.waterflows_subglacial(surfdem, beddem, dx; gamma=WWFS.GAMMA)

length(out.routing.sinks), sum(out.pressmelt.sc_locs)
```

Set `gamma=0` to disable deflection/supercooling effects and recover routing
equivalent to potential-only behavior.

## Output structure

`waterflows_subglacial` returns a nested named tuple with sections:

- `routing`: routed fields (`area`, `dir`, `sinks`, `pits`, ...)
- `pressmelt`: supercooling flags and deflection diagnostics
- `lakes`: fixed/free-surface lake-depth diagnostics
- `sink_catchments`: optional masks and per-sink fluxes

## Key controls

- `gamma`: supercooling/deflection strength (`0` disables this effect)
- `avoid_sc`: whether supercooled cells act as barriers
- `floatfrac`: flotation fraction (scalar or field)
- `source`: source term to route
- `mask`: active routing mask

Other useful controls include `ctch_sinks` for grouped sink-flux diagnostics,
and `drain_pits`, `bnd_as_sink`, `nan_as_sink` inherited from core routing.

## Suggested examples

- `examples/wwfs-simple.jl`: quick intro
- `examples/subglacially/valley-glacier.jl`: valley setup
- `examples/subglacially/ice-sheet-margin-shmip.jl`: SHMIP-like geometry

You can run those from the examples environment, for example:

```julia
include("subglacially/valley-glacier.jl")
```

## Physical notes

- Inputs are expected on a regular Cartesian grid with scalar `dx`.
- Routing assumes local D8 connectivity.
- Some diagnostics (e.g. supercooling masks) are model approximations and
  should be interpreted accordingly.

See also: [Tutorial](@ref), [Examples](@ref), and [API Reference](@ref).
