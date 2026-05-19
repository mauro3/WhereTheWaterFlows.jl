# [Tutorial](@id Tutorial)

This tutorial walks through a complete flow-routing analysis on a synthetic DEM.

## Setup

```julia
using WhereTheWaterFlows, GLMakie   # GLMakie (or CairoMakie) enables the plt_* functions
const WWF = WhereTheWaterFlows
```

## Build a synthetic DEM

We use a small function that produces a landscape with a few hills, valleys, and
a deliberate closed depression (pit):

```julia
function peaks2(n=100, randfac=0.05)
    coords = range(-π, π, length=n)
    return coords, coords,
           sin.(coords) .* cos.(coords') .-
           0.7 .* (sin.(coords .+ 1) .* cos.(coords')).^8 .+
           randfac .* randn(n, n)
end

x, y, dem = peaks2(200)
```

## Run the flow router

```julia
area, slen, dir, nout, nin, sinks, pits, c, bnds = waterflows(dem)
```

`waterflows` returns ten values as a `NamedTuple`:

| Name | Description |
|---|---|
| `area` | Upslope area at each cell (number of cells draining into it, or a physical flux if `cellarea` is supplied) |
| `slen` | Stream length: number of cells from the farthest upstream source |
| `dir` | D8 flow direction at each cell, encoded as an integer 1–9 |
| `nout` | 1 if the cell has a downstream neighbour, 0 if it is a pit |
| `nin` | Number of upstream neighbours flowing into each cell |
| `sinks` | `Vector{CartesianIndex{2}}` of boundary/NaN-adjacent sink cells |
| `pits` | `Vector{CartesianIndex{2}}` of interior pit cells (local minima) |
| `c` | Integer catchment map; each pit or sink has its own colour |
| `bnds` | `Vector` of boundary-cell lists, one per catchment |
| `flowdir_extra_output` | Extra output from a custom `flowdir_fn` (nothing for the default) |

## Visualise with Makie

```julia
# log10 upslope area with pit locations marked
plt_area(x, y, area, pits)

# colour-coded catchment map
plt_catchments(x, y, c)

# stream length
heatmap(x, y, slen)
```

## Delineate a single catchment

Pick the cell with the largest upslope area along row 50:

```julia
i = 50
j = findmax(area[i, :])[2]

cc = catchment(dir, CartesianIndex(i, j))   # BitArray, true inside the catchment
heatmap(x, y, cc)
scatter!(x[i], y[j], markersize=20)
```

## Fill depressions

After calling `waterflows`, pits have already been *routed through* (by default
`drain_pits=true`), so the flow directions are consistent across depressions.
`fill_dem` raises every pit cell to the elevation of its spillway, useful for
visualising lake depth or for further analysis:

```julia
demf = fill_dem(dem, sinks, dir)

# lake depth
heatmap(x, y, demf .- dem)
```

## Key options of `waterflows`

| Keyword | Default | Meaning |
|---|---|---|
| `drain_pits` | `true` | Route water through pits by reversing flow to the lowest spillway |
| `bnd_as_sink` | `false` | Treat domain boundary cells as sinks |
| `nan_as_sink` | `true` | Cells adjacent to NaN cells become sinks |
| `cellarea` | ones array | Source flux per cell; can be a tuple for multiple tracers |
| `feedback_fn` | `nothing` | Applied to accumulated area at each cell before routing further downstream |

## Routing physical fluxes

Set `cellarea` to a physical value (e.g. precipitation in m³/s per cell) to
accumulate real fluxes rather than cell counts:

```julia
precip = 1e-3 .* ones(size(dem))          # uniform 1 mm/s
discharge, = waterflows(dem, precip)            # area now in m³/s
```

Pass a tuple to route several quantities simultaneously:

```julia
water  = ones(size(dem))
tracer = zeros(size(dem)...)
tracer[[40, 678, 4560, 13476]] .= 1
(water_area, tracer_area), slen, dir, = waterflows(dem, (water, tracer))
```

## Subglacial routing example

The `examples/` directory contains two more advanced scripts:

- `examples/subglacial-routing.jl` — route subglacial water under an ice sheet
  using a hydraulic potential instead of a surface DEM.
- `examples/subglacial-routing-feedback.jl` — the same but with a `feedback_fn`
  that limits flux to non-negative values, demonstrating self-feedback.
