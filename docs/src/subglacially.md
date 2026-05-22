# [Subglacially](@id SubglaciallyGuide)

`WhereTheWaterFlows.Subglacially` extends core routing for subglacial settings.

Compared with plain `waterflows`, it provides:

- routing based on Shreve hydraulic potential,
- pressure-melting/supercooling-aware flow deflection,
- lake-depth diagnostics,
- optional per-sink flux aggregation.

This module is aimed at routing water under ice using the Shreve hydraulic
potential.  Optionally, pressure-melting-point effects (the Röthlisberger
deflection) are accounted for.

## Hydraulic potential

The Shreve hydraulic potential φ used for routing is

```
φ = f · H · (ρᵢ/ρ_w) + (z_s − H)
```

where *H* is ice thickness, *z_s* is surface elevation, *f* is the flotation
fraction, and ρᵢ, ρ_w are ice and water density.  At full flotation (*f* = 1)
this is the standard Shreve potential.  The bed elevation is z_b = z_s − H.

Water flows down the gradient of φ, not down the gradient of the bed.

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

Set `gamma=0` to disable deflection/supercooling effects and recover
potential-only behaviour.

## Output structure

`waterflows_subglacial` returns a nested named tuple with four top-level keys:

### `routing`

All fields from the core `waterflows` call, plus `phi` (the hydraulic
potential used for routing):

- `area`: named tuple with four routed fields:
  - `total`: total accumulated water flux (source + melt; in m³/s when `source` is in m/s)
  - `extra`: accumulated flux from melt only (dissipation + pressure melting)
  - `dissipation_melt_rate`: local melt rate from energy dissipation [m/s]
  - `pressure_melt_rate`: local melt rate from pressure-melting-point effects [m/s]
- `slen, dir, nout, nin, sinks, pits, c, bnds`: same as `waterflows` (see [Tutorial](@ref))
- `phi`: the Shreve hydraulic potential field

### `pressmelt`

Diagnostics from the flow-deflection step:

- `sc_locs`: Boolean mask of cells where supercooling occurs (outflow edge is supercooled)
- `kappas`: deflection angle in units of π/4 at each cell
- `dir_og`: flow directions before deflection

### `lakes`

Lake-depth diagnostics derived from `fill_dem` applied to the hydraulic potential:

- `depth_fixed_surface`: lake depth assuming the ice surface is fixed [m water equivalent]
- `depth_free_surface`: lake depth assuming the ice surface adjusts freely [m]

### `sink_catchments`

Only populated when `ctch_sinks` is non-empty:

- `masks`: vector of `BitMatrix`, one per sink set
- `fluxes`: named tuple `(total, dissipation, pressmelt)`, each a vector of
  total flux values (m³/s) for the corresponding sink set

## Key controls

| Argument / Keyword | Default | Meaning |
|--------------------|---------|---------|
| `surfdem`, `beddem` | — | Surface and bed elevation arrays (same size) |
| `dx` | — | Grid spacing in metres (must be equal in x and y) |
| `floatfrac` | `1` | Flotation fraction (scalar or array); `1` = full flotation |
| `source` | `ones(size(surfdem))` | Meltwater input per unit area [m/s]; multiply by `dx²` to get m³/s per cell |
| `mask` | all `true` | Active routing mask; `false` cells are set to NaN in φ |
| `gamma` | `GAMMA` (≈−0.31) | Röthlisberger constant controlling deflection strength; set to `0` to disable |
| `avoid_sc` | `false` | If `true`, supercooled cells become barriers (mass is lost there) |
| `ctch_sinks` | `[]` | Sink-area sets for per-catchment flux diagnostics |
| `rhow`, `rhoi` | 1000, 910 | Water and ice densities [kg m⁻³] |
| `drain_pits`, `bnd_as_sink`, `nan_as_sink` | `true` | Inherited from core routing |

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
- Supercooling diagnostics (`sc_locs`) are a model approximation and should be
  interpreted with care.
- Negative `dissipation_melt_rate` can occur where the breach algorithm routes
  water slightly uphill out of a filled depression; the corresponding freeze
  upstream balances this.

See also: [Tutorial](@ref), [Examples](@ref), and [API Reference](@ref).
