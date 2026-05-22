# [Randomly](@id RandomlyGuide)

`WhereTheWaterFlows.Randomly` provides Monte Carlo helpers for uncertainty
propagation around deterministic routing models.

The main workflow is:

1. define uncertainty for uncertain fields,
2. construct `model`, `sample`, and `reduce!` with `make_fns_*`,
3. run many realizations with `map_mc`.

## Core pieces

- `Uncertainty`: uncertainty model for a field
- `make_fns_subaerial`: stochastic wrapper around `WhereTheWaterFlows.waterflows`
- `make_fns_subglacial`: stochastic wrapper around
  `WhereTheWaterFlows.Subglacially.waterflows_subglacial`
- `map_mc`: Monte Carlo loop with aggregation

## Minimal subaerial run

```@example randomly
using WhereTheWaterFlows
using Random, Statistics

Random.seed!(42)

WWFR = WhereTheWaterFlows.Randomly

n = 80
dx = 100.0
x = range(-pi, pi, length=n)
dem = sin.(x) .* cos.(x')
source = fill(1e-3 / dx^2, size(dem))

dem_uc = WWFR.Uncertainty()  # keep DEM fixed
source_uc = WWFR.Uncertainty(absuc=0.0, reluc=0.2, correlation_length=15 * dx)

ctch_sinks = [CartesianIndices((2:2, 2:n-1))[:]]

model, sample, reduce! = WWFR.make_fns_subaerial(dx,
                                                 dem, dem_uc,
                                                 source, source_uc,
                                                 ctch_sinks)

aggr = WWFR.map_mc(model, sample, reduce!, 6; progressmeter=false)

mean(aggr.catchment_fluxes[1]), std(aggr.catchment_fluxes[1])
```

## Interpreting outputs

For subaerial runs, the aggregate contains e.g.:

- `areas_total`: mean routed area/discharge field
- `stream_length`: mean stream-length field
- `catchments`: per-sink occurrence frequency (0-1)
- `catchment_fluxes`: one vector of sample fluxes per sink set

For subglacial runs, additional fields are included (lakes, supercooling,
pressure-melt diagnostics).

## Suggested examples

- `examples/wwfr-simple.jl`: quick intro
- `examples/randomly/source-uncertainty-sweep.jl`: sensitivity-style sweep

See also: [Examples](@ref) and [API Reference](@ref).
