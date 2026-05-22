# [Subglacially](@id SubglaciallyGuide)

`WhereTheWaterFlows.Subglacially` extends core routing for subglacial settings.

Compared with plain `waterflows`, it provides:

- routing based on hydraulic potential,
- pressure-melting/supercooling-aware deflection,
- lake-depth diagnostics,
- optional sink-catchment flux aggregation.

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

## Suggested examples

- `examples/wwfs-simple.jl`: quick intro
- `examples/subglacially/valley-glacier.jl`: valley setup
- `examples/subglacially/ice-sheet-margin-shmip.jl`: SHMIP-like geometry

See also: [Tutorial](@ref), [Examples](@ref), and [API Reference](@ref).
