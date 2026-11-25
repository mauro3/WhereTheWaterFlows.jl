# [WhereTheWaterFlows.Randomly (WWFR)](@id RandomlyGuide)

`WWFR` provides Monte Carlo wrappers to propagate field
uncertainty through deterministic routing models.

It supports subaerial routing uncertainty (`make_fns_subaerial`),
subglacial routing uncertainty (`make_fns_subglacial`). It implements the spatial
uncertainty in the input fields using configurable Gaussian random fields (GRFs)
which produce spatially correlated noise.

## Workflow

Typical usage is:

1. define uncertainty models for uncertain fields,
2. build `model`, `sample`, `initialize`, `aggregate!`, and `finalize!` with `make_fns_*`,
3. run Monte Carlo with `map_mc`.

The key design idea is that `map_mc` is generic: it does not know anything
about routing fields. It calls five functions you provide:

- `sample()`: generate one stochastic input realization,
- `model(args...)`: run one deterministic forward model,
- `initialize()`: create fresh aggregate storage,
- `aggregate!(aggr, output)`: accumulate one realization in place,
- `finalize!(aggr)`: normalize aggregate storage in place.

Both wrappers return the named tuple `(; model, sample, initialize, aggregate!, finalize!)`.
Use named destructuring to select the callbacks needed by your workflow.

## Uncertainty model

`Uncertainty` stores how one field should be perturbed:

```julia
Uncertainty(; absuc=0,
              reluc=0,
              correlation_length=1,
              covariance_fn=gaussian_kernel,
              abs_bounds=(-Inf, Inf))
```

For a base field `f`, each realization is generated as:

1. draw a zero-mean unit-variance correlated GRF `xi`,
2. scale it with absolute and relative amplitudes,
3. optionally clip perturbation values,
4. add perturbation back to the field.

In code terms (as implemented):

```julia
delta = xi .* absuc .+ xi .* f .* reluc
delta = clamp.(delta, abs_bounds...)
f_realization = f .+ delta
```

### Parameter meaning

| Field | Meaning |
|---|---|
| `absuc` | Absolute perturbation amplitude (scalar or array-like) |
| `reluc` | Relative perturbation amplitude, scaled by local field value |
| `correlation_length` | Correlation length in physical units (same units as `dx`) |
| `covariance_fn` | Spatial kernel (`gaussian_kernel`, `exponential_kernel`, or custom) |
| `abs_bounds` | Bounds applied to perturbation `delta` before adding to field |

## GRFs: how they work

WWFR uses FFT-based GRF sampling (following Raess et al., 2019).

Correlation length conversion: if the grid spacing is `dx`, then WWFR converts:
```julia
len = correlation_length / dx
```
so kernels operate in cell units.

Kernel choice:
- `gaussian_kernel`: smoother, more diffuse perturbations
- `exponential_kernel`: rougher, shorter-scale structure for same nominal length
- custom kernel: pass any `kernel_fn(nx, ny, len)` compatible function

## [Seeded sampling and replay](@id SeededSampling)

By default, samples are drawn from Julia's global RNG, so a single draw can only be
reproduced by replaying the whole sequence. Alternatively, pass a table of per-draw,
per-field seeds to `map_mc`; any draw can then be regenerated on its own.

```@example seeded_sampling
using WhereTheWaterFlows, Random
const WWF = WhereTheWaterFlows
const WWFR = WWF.Randomly

nx, ny, dx = 24, 22, 100.0
dem = [1000.0 - 5i - j for i in 1:nx, j in 1:ny]
source = fill(1e-8, nx, ny)
dem_uc = WWFR.Uncertainty()
source_uc = WWFR.Uncertainty(reluc=0.1, correlation_length=2dx)
ctch_sinks = [WWF.waterflows(dem).sinks]

(; model, sample, initialize, aggregate!, finalize!) = WWFR.make_fns_subaerial(
    dx, dem, dem_uc, source, source_uc, ctch_sinks)
seeds = WWFR.make_sample_seeds(8, (:dem, :source); seed=42)
aggr = WWFR.map_mc(model, sample, initialize, aggregate!, finalize!, 8;
                   seeds, progressmeter=false)

# Replay draw 3 on its own
draw = sample(; seeds=aggr.sample_seeds[3])
replayed = model(draw.fields...)

(; ids=aggr.sample_ids, rng_type=aggr.rng_type)
```

In seeded mode `sample(; seeds)` returns a
[`WWFR.SeededSample`](@ref WhereTheWaterFlows.Randomly.SeededSample) holding
`.fields` (the model arguments), `.seeds` and `.rng_type`. Each field gets its own
`Random.Xoshiro`. The seed-record keys of the built-in wrappers are:

| Wrapper | Seed record keys |
|---|---|
| `make_fns_subaerial` | `(:dem, :source)` |
| `make_fns_subglacial` | `(:surface, :bed, :floatfrac, :source)` |

`initialize()` must then return a named tuple; `map_mc` returns it extended with
`sample_ids`, `sample_seeds` (both in accumulation order, matching the per-sample flux
vectors) and `rng_type`. These are not touched by `finalize!`.

To run a subset or a reordered batch, pass the selected seeds and their original IDs:

```@example seeded_sampling
all_seeds = WWFR.make_sample_seeds(1000, (:dem, :source); seed=42)
ids = [73, 1000, 12]
selected = WWFR.map_mc(model, sample, initialize, aggregate!, finalize!, length(ids);
    seeds=all_seeds[ids], sample_ids=ids, progressmeter=false)
selected.sample_ids
```

A field's seeds depend only on the master seed and the field name. So changing `n`, or
adding, removing or reordering other fields, keeps them; this allows paired comparisons
of runs that differ in one uncertain field. The same seed gives the same field only for
the same grid, kernel, correlation length, precision and code/package versions, so
record those along with the seeds. Zero uncertainty consumes no random numbers.

A custom sampler opts in by accepting `sample(; seeds)` and returning a
`SeededSample`; the model, `aggregate!` and `finalize!` are unchanged:

```@example custom_seeded_sampling
using WhereTheWaterFlows, Random
const WWF = WhereTheWaterFlows
const WWFR = WWF.Randomly

function sample(; seeds)
    rng = MersenneTwister(seeds.forcing)
    return WWFR.SeededSample((randn(rng),), seeds; rng_type=MersenneTwister)
end
initialize() = (; values=Float64[])
aggregate!(aggr, value) = push!(aggr.values, value)
finalize!(aggr) = nothing
seeds = WWFR.make_sample_seeds(3, (:forcing,); seed=42)
aggr = WWFR.map_mc(identity, sample, initialize, aggregate!, finalize!, 3;
                   seeds, progressmeter=false)
@assert aggr.values[2] == randn(aggr.rng_type(aggr.sample_seeds[2].forcing))
```

## What can be configured

High-level Monte Carlo controls:
- `map_mc(model, sample, initialize, aggregate!, finalize!, n;
  progressmeter=true, seeds=nothing, sample_ids=1:n)`
- `n`: number of realizations (`1 <= n <= 2048` enforced)
- `progressmeter`: set `false` for quiet batch/CI runs
- `seeds`: one per-field seed record per draw, enabling seeded mode
- `sample_ids`: one integer ID per draw, used for seeded metadata (not RNG derivation)

Subaerial wrapper (`make_fns_subaerial`):
- uncertain fields: `dem`, `source`
- routing controls: `drain_pits`, `bnd_as_sink`, `nan_as_sink`
- sink groups for aggregation: `ctch_sinks`
- catchment storage: `catchment_storage=:blocksparse` (default), or `:dense`
- catchment tile size: `catchment_blocksize=50` (square tiles, `:blocksparse` only)

Subglacial wrapper (`make_fns_subglacial`):
- uncertain fields: `surfdem`, `beddem`, `floatfrac`, `source`
  - Note: unlike `waterflows_subglacial`, which accepts scalar `floatfrac`,
    `make_fns_subglacial` requires `floatfrac` to be an array.
- physical controls: `gamma`, `rhow`, `rhoi`
- processing controls: `mask`, `ctch_sinks`, `min_lake_depth`
- catchment storage: `catchment_storage=:blocksparse` (default), or `:dense`
- catchment tile size: `catchment_blocksize=50` (square tiles, `:blocksparse` only)

## Why reduction is needed

A Monte Carlo run can easily produce hundreds of large 2D outputs. Storing all
realizations is expensive and often unnecessary. Aggregation keeps only summary
statistics (means, frequencies, selected per-sample vectors), which is both
memory-efficient and directly useful for interpretation.

In WWFR wrappers, reduction happens in three stages:

1. `initialize()` initializes aggregate storage.
2. `aggregate!(aggr, sample_output)` updates aggregate values per realization.
3. `finalize!(aggr)` finalizes (typically divides accumulated maps by `n`).

`map_mc` calls initialization and finalization once each and returns the aggregate
storage (in seeded mode extended by the seed metadata). Return values from the two
mutating callbacks are ignored.

## How outputs are reduced

### Catchment storage

Both wrappers return `aggr.catchments` as a three-dimensional Float16 array with
shape `(nx, ny, n_groups)`. The default is a `BlockArrays.BlockArray` containing
50×50×1 tiles (smaller at spatial edges). Untouched tiles share lazy zeros;
only tiles reached by a catchment are allocated. Thus for cases of catchments which 
occupy only a small-ish fraction of the domain, this sparse storage can decrease
 memory requirements considerably. Counts are accumulated in this
storage and normalized in place, so final probabilities remain tiled. Tile support
is the union over all samples, not just the last realization.

Set `catchment_blocksize` in either `make_fns_*` call to use square tiles of a
different size, for example `catchment_blocksize=25`. The default remains 50.
The block size only affects memory use and update performance, not the numerical
results. It is ignored when `catchment_storage=:dense`.

Use ordinary indexing, such as `aggr.catchments[x, y, group]`, and
`view(aggr.catchments, :, :, group)` for a zero-copy matrix view. Implicit zero is
valid zero membership, not missing data. Groups remain independent and may overlap.
The deterministic `Subglacially` output still contains the same packed Boolean masks.

Set `catchment_storage=:dense` in either `make_fns_*` call to obtain a conventional
`Array{Float16,3}` instead. This changes only membership storage, not sampling,
routing, flux records, other aggregate fields, or precision. Float16 counts are exact
through 2048 samples; final probabilities are rounded to Float16 as before.

The tiled output is intended for reading. A nonzero scalar write into an untouched
lazy tile is not supported; use `Array(aggr.catchments)` if a conventional editable
array is needed. Slicing, generic array operations and plotting can materialize dense
storage: crop/stride views before converting or plotting large grids. The reducer's
updates and finalization explicitly avoid these materializations.

### Subaerial

The `aggregate!` function returned by `make_fns_subaerial` accumulates:
- `areas_total += output.area`
- `stream_length += output.slen`
- `catchments[:, :, i] += catchment(output.dir, ctch_sinks[i])`
- `catchment_fluxes[i]` stores one scalar per sample (not averaged in place)

At the finalize step map-style fields are normalized by `n_samples` (so they are
means/frequencies) and flux vectors remain per-sample records by design.

### Subglacial

The `aggregate!` function returned by `make_fns_subglacial` accumulates:
routed maps and diagnostics (`areas_total`,
`lake_depth_*`, `lake_mask_*`, `sc_locs`, `kappas`, `catchments`) and appends
sink-flux records to vectors (`catchment_fluxes.total/dissipation/pressmelt`).

Catchments may overlap. Membership is accumulated only at true mask cells, using
Float16 counts and in-place normalization. Each group retains a full-grid BitMatrix
per realization; the aggregate uses the catchment storage described above.

At the finalize step map-style fields are normalized by `n_samples` (so they are
means/frequencies) and flux vectors remain per-sample records by design.

When `ctch_sinks` is non-empty, `aggr.catchment_fluxes` has subglacial-specific
structure:
- `aggr.catchment_fluxes` is a named tuple
  `(total, dissipation, pressmelt)`.
- Each field is a `Vector` with length `length(ctch_sinks)`.
- Element `i` of each field is a `Vector{Float32}` with length `n_samples`,
  containing the per-realization flux for outlet group `i`.

So for subglacial runs, the total-flux distribution at outlet group `i` is
accessed as `aggr.catchment_fluxes.total[i]`.
This differs from the subaerial wrapper, where per-outlet total fluxes are
stored directly as `aggr.catchment_fluxes[i]` (no `.total`).

Minimal usage pattern:

```julia
using Statistics

for i in eachindex(ctch_sinks)
    fluxes = aggr.catchment_fluxes.total[i]   # Vector{Float32}, length n_samples
    println("Outlet $i: mean=$(mean(fluxes)), std=$(std(fluxes))")
end
```

## `ctch_sinks`: what it means

`ctch_sinks` defines *sink groups*, as described in [Subglacially: Defining outlet groups](@ref "Defining outlet groups (`ctch_sinks`)"), for which WWFR reports catchment membership
and integrated fluxes.

Each element of `ctch_sinks` is a collection of sink cells (typically a
`Vector{CartesianIndex}` or `CartesianIndices`). For each Monte Carlo sample,
WWFR computes the full upstream catchment draining to each sink group.
You can pass multiple groups, e.g. upper/lower terminus sectors, to compare
how uncertainty redistributes flux between outlets.

This enables statistics like:

- probability that a cell drains to outlet group `i` (`aggr.catchments[:, :, i]`),
- distribution of total flux entering outlet group `i`
  (`aggr.catchment_fluxes[i]` in subaerial, and
  `aggr.catchment_fluxes.total[i]` in subglacial).


## Defining custom `make_fns_*`

The built-in `make_fns_subaerial`/`make_fns_subglacial` are templates. If your
model has different outputs, constraints, or diagnostics, define your own callbacks:

- `model(args...)`
- `sample()`
- `initialize()`
- `aggregate!(aggr, out)`
- `finalize!(aggr)`

then call `map_mc(model, sample, initialize, aggregate!, finalize!, n)`.

Minimal skeleton:

```julia
model(a, b) = my_forward_model(a, b)
sample() = (draw_a(), draw_b())

function initialize()
    return (sumfield = Ref(0.0), n = Ref(0))
end

function aggregate!(aggr, out)
    aggr.sumfield[] += out.metric
    aggr.n[] += 1
    return aggr
end

function finalize!(aggr)
    aggr.sumfield[] /= aggr.n[]
    return aggr
end
```

`initialize` creates storage once at the beginning of the Monte Carlo run.
`aggregate!` updates it after each forward model evaluation, and `finalize!`
normalizes the result at the end. Use mutable storage (arrays, `Ref`s, or a mutable
struct) for values updated in place.

`map_mc` currently runs samples sequentially.

## Minimal subaerial example

```@example randomly
using WhereTheWaterFlows, CairoMakie
using Statistics
using Random; Random.seed!(42)

WWFR = WhereTheWaterFlows.Randomly

n = 80
dx = 100.0
x = range(-pi, pi, length=n)
dem = sin.(x) .* cos.(x')
source = fill(1e-3 / dx^2, size(dem))

dem_uc = WWFR.Uncertainty()  # fixed DEM
source_uc = WWFR.Uncertainty(absuc=0.0, reluc=0.2, correlation_length=15 * dx)

ctch_sinks = [CartesianIndices((1:n, 1:1))[:]]

(; model, sample, initialize, aggregate!, finalize!) = WWFR.make_fns_subaerial(dx,
                                                 dem, dem_uc,
                                                 source, source_uc,
                                                 ctch_sinks)

aggr = WWFR.map_mc(model, sample, initialize, aggregate!, finalize!, 100; progressmeter=false)

mu, sd = mean(aggr.catchment_fluxes[1]), std(aggr.catchment_fluxes[1])
hist(aggr.catchment_fluxes[1],
     axis=(title  = "mean = $(round(mu, sigdigits=4)), std = $(round(sd, sigdigits=4))",
     xlabel = "flux", ylabel="freq."))
```

## Output interpretation

### Subareal aggregate outputs (`make_fns_subaerial`)
For subaerial runs, aggregated outputs:

- `areas_total`: Monte Carlo mean routed area/discharge field
- `stream_length`: Monte Carlo mean stream-length field
- `catchments`: per-sink occurrence frequency map (0-1)
- `catchment_fluxes`: per-sink vectors of sample fluxes

### [Subglacial aggregate outputs (`make_fns_subglacial`)](@id SubglAggr)

The table below consolidates all fields produced by `aggr = map_mc(...)` when
using `make_fns_subglacial`.

| Field | Shape / type | Meaning | Reduction semantics |
|---|---|---|---|
| `aggr.areas_total` | `Matrix{Float32}` | Total routed discharge proxy [m^3/s] | Mean over samples |
| `aggr.areas_extra` | `Matrix{Float32}` | Extra routed discharge from subglacial melt terms [m^3/s] | Mean over samples |
| `aggr.melt_rate` | `Matrix{Float32}` | Dissipation + pressure-melt generation rate [m/s] | Mean over samples |
| `aggr.lake_depth_fixed_surface` | `Matrix{Float32}` | Lake depth for fixed-surface assumption [m] | Mean over samples |
| `aggr.lake_mask_fixed_surface` | `Matrix{Float32}` | Lake-occurrence frequency (depth `> min_lake_depth`) for fixed-surface case [0, 1] | Fraction over samples |
| `aggr.lake_depth_free_surface` | `Matrix{Float32}` | Lake depth for free-surface assumption [m] | Mean over samples |
| `aggr.lake_mask_free_surface` | `Matrix{Float32}` | Lake-occurrence frequency (depth `> min_lake_depth`) for free-surface case [0, 1] | Fraction over samples |
| `aggr.sc_locs` | `Matrix{Float32}` | Supercooling occurrence frequency [0, 1] | Fraction over samples |
| `aggr.kappas` | `Matrix{Float32}` | Mean supercooling deflection | Mean over samples |
| `aggr.catchments` | `BlockArray{Float16,3}` by default (`nx x ny x n_groups`, `catchment_blocksize`-sized tiles); `Array{Float16,3}` with `catchment_storage=:dense` | Catchment-membership frequency for each `ctch_sinks` group [0, 1] | Fraction over samples |
| `aggr.catchment_fluxes.total` | `Vector{Vector{Float32}}` (length `n_groups`) | Per-group total outlet flux [m^3/s] | Per-sample vectors (not averaged) |
| `aggr.catchment_fluxes.dissipation` | `Vector{Vector{Float32}}` (length `n_groups`) | Per-group dissipation-melt flux contribution records [m^3/s] | Per-sample vectors (not averaged) |
| `aggr.catchment_fluxes.pressmelt` | `Vector{Vector{Float32}}` (length `n_groups`) | Per-group pressure-melt flux contribution records [m^3/s] | Per-sample vectors (not averaged) |
| `aggr.lake_vol_fixed_surface` | `Vector{Float32}` (length `n_samples`) | Lake volume [m^3]: `dx^2` times the sum of full depths `> min_lake_depth` (fixed-surface case) | Per-sample vector (not averaged) |
| `aggr.lake_vol_free_surface` | `Vector{Float32}` (length `n_samples`) | Lake volume [m^3]: `dx^2` times the sum of full depths `> min_lake_depth` (free-surface case) | Per-sample vector (not averaged) |
| `aggr.n_samples[]` | `Int` (`Ref`) | Number of Monte Carlo realizations accumulated | Final count |

Notes:

- Seeded runs also contain `sample_ids`, `sample_seeds` and `rng_type`, see
  [Seeded sampling and replay](@ref SeededSampling).
- For `aggr.catchment_fluxes.*`, element `i` is the flux distribution for outlet
  group `i`; each inner vector has length `aggr.n_samples`.
- If `ctch_sinks` is empty, `aggr.catchments` has zero groups and
  `aggr.catchment_fluxes.*` are empty vectors.

## Practical guidance

- Use `absuc` for additive error floors; `reluc` for proportional uncertainty.
- Keep `correlation_length` meaningfully smaller than domain size (as a rule of
  thumb, at least around 3x smaller than domain extent).
- Use [seed records](@ref SeededSampling) for reproducible runs and replay of single
  draws.
- Start with small `n` while validating model setup; scale up after checks.

## Examples

- `examples/wwfr-simple.jl`: quick start
- `examples/randomly/source-uncertainty-sweep.jl`: source-vs-DEM uncertainty
  comparison and sensitivity sweep
- `examples/subglacially/ice-cap-full-workflow.jl`: subglacial deterministic+
  Monte Carlo workflow, including sink discovery and `ctch_sinks`
  partitioning for per-outlet flux distributions

From the `examples/` environment:

```julia
include("wwfr-simple.jl")
include("randomly/source-uncertainty-sweep.jl")
include("subglacially/ice-cap-full-workflow.jl")
```

See also: [Examples](@ref ExamplesPage), [Subglacially](@ref SubglaciallyGuide), and
[API Reference](@ref APIReference).
