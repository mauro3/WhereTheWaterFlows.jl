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
2. build `model`, `sample`, and `reduce!` with `make_fns_*`,
3. run Monte Carlo with `map_mc`.

The key design idea is that `map_mc` is generic: it does not know anything
about routing fields. It only calls three functions you provide:

- `sample()`: generate one stochastic input realization,
- `model(args...)`: run one deterministic forward model,
- `reduce!`: aggregate outputs over many realizations.

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

WWFR uses FFT-based GRF sampling (following Raess et al., 2019):

1. sample white noise `Z`,
2. transform with FFT,
3. multiply by `sqrt(FFT(kernel))` to enforce target covariance,
4. inverse FFT to get a correlated field.

Internally, the domain is padded before FFT to improve non-periodic behavior
and to reach FFT-friendly sizes (factors of 2, 3, 5, 7).

### Correlation length conversion

If the grid spacing is `dx`, then WWFR converts:

```julia
len = correlation_length / dx
```

so kernels operate in cell units.

### Kernel choice

- `gaussian_kernel`: smoother, more diffuse perturbations
- `exponential_kernel`: rougher, shorter-scale structure for same nominal length
- custom kernel: pass any `kernel_fn(nx, ny, len)` compatible function

## What can be configured

### High-level Monte Carlo controls

- `map_mc(model, sample, reduce!, n; progressmeter=true)`
- `n`: number of realizations (`n <= 2048` enforced)
- `progressmeter`: set `false` for quiet batch/CI runs

### Subaerial wrapper (`make_fns_subaerial`)

- uncertain fields: `dem`, `source`
- routing controls: `drain_pits`, `bnd_as_sink`, `nan_as_sink`
- sink groups for aggregation: `ctch_sinks`

### Subglacial wrapper (`make_fns_subglacial`)

- uncertain fields: `surfdem`, `beddem`, `floatfrac`, `source`
- physical controls: `gamma`, `rhow`, `rhoi`
- processing controls: `mask`, `ctch_sinks`, `min_lake_depth`

## Why reduction is needed

A Monte Carlo run can easily produce hundreds of large 2D outputs. Storing all
realizations is expensive and often unnecessary. `reduce!` keeps only summary
statistics (means, frequencies, selected per-sample vectors), which is both
memory-efficient and directly useful for interpretation.

In WWFR wrappers, reduction happens in three stages:

1. `reduce!()` initializes aggregate storage.
2. `reduce!(aggr, sample_output)` updates aggregate values per realization.
3. `reduce!(aggr)` finalizes (typically divides accumulated maps by `n`).

## How outputs are reduced

### `make_fns_subaerial`

Per sample, it accumulates:

- `areas_total += output.area / dx^2`
- `stream_length += output.slen`
- `catchments[:, :, i] += catchment(output.dir, ctch_sinks[i])`
- `catchment_fluxes[i]` stores one scalar per sample (not averaged in place)

At finalize step it divides `areas_total`, `stream_length`, and `catchments` by
`n_samples`, so these become Monte Carlo means/frequencies.

### `make_fns_subglacial`

Per sample, it accumulates routed maps and diagnostics (`areas_total`,
`lake_depth_*`, `lake_mask_*`, `sc_locs`, `kappas`, `catchments`) and appends
sink-flux records to vectors (`catchment_fluxes.total/dissipation/pressmelt`).

At finalize step it normalizes map-style fields by `n_samples` (so they are
means/frequencies). Flux vectors remain per-sample records by design.

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

`ctch_sinks` defines *sink groups* for which WWFR reports catchment membership
and integrated fluxes.

Each element of `ctch_sinks` is a collection of sink cells (typically a
`Vector{CartesianIndex}` or `CartesianIndices`). For each Monte Carlo sample,
WWFR computes the full upstream catchment draining to each sink group.

This enables statistics like:

- probability that a cell drains to outlet group `i` (`aggr.catchments[:, :, i]`),
- distribution of total flux entering outlet group `i`
  (`aggr.catchment_fluxes[i]` in subaerial, and
  `aggr.catchment_fluxes.total[i]` in subglacial).

Example (left-margin outlet band):

```@example randomly
n = 80
ctch_sinks = [CartesianIndices((2:2, 2:n-1))[:]]
length(ctch_sinks), length(ctch_sinks[1])
```

You can pass multiple groups, e.g. upper/lower terminus sectors, to compare
how uncertainty redistributes flux between outlets.

## Defining custom `make_fns_*`

The built-in `make_fns_subaerial`/`make_fns_subglacial` are templates. If your
model has different outputs, constraints, or diagnostics, define your own trio:

- `model(args...)`
- `sample()`
- `reduce!` (three-method interface)

then call `map_mc(model, sample, reduce!, n)`.

Minimal skeleton:

```julia
model(a, b) = my_forward_model(a, b)
sample() = (draw_a(), draw_b())

function reduce!()
    return (sumfield = 0.0, n = Ref(0))
end

function reduce!(aggr, out)
    aggr.sumfield += out.metric
    aggr.n[] += 1
    return aggr
end

function reduce!(aggr)
    aggr.sumfield /= aggr.n[]
    return aggr
end
```

The `reduce!` function has three methods: 0-arg method sets up the storage needed
in the reduction (this gets called once at the beginning of the mc-iterations); the
2-arg method is called after each forward model evaluation and reduces/aggregates the
forward model output into what is stored; the 1-arg method is then called at the end
to finalize the aggregated results.

Note: only `reduce!` is allowed to not be thread-safe, `model` and `sample` need to be thread-safe.

This pattern is often the cleanest route for problem-specific uncertainty
metrics.

## Minimal subaerial example

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

dem_uc = WWFR.Uncertainty()  # fixed DEM
source_uc = WWFR.Uncertainty(absuc=0.0, reluc=0.2, correlation_length=15 * dx)

ctch_sinks = [CartesianIndices((2:2, 2:n-1))[:]]

model, sample, reduce! = WWFR.make_fns_subaerial(dx,
                                                 dem, dem_uc,
                                                 source, source_uc,
                                                 ctch_sinks)

aggr = WWFR.map_mc(model, sample, reduce!, 6; progressmeter=false)

mean(aggr.catchment_fluxes[1]), std(aggr.catchment_fluxes[1])
```

## Output interpretation

For subaerial runs, aggregated outputs include:

- `areas_total`: Monte Carlo mean routed area/discharge field
- `stream_length`: Monte Carlo mean stream-length field
- `catchments`: per-sink occurrence frequency map (0-1)
- `catchment_fluxes`: per-sink vectors of sample fluxes

For subglacial runs, aggregation also includes lake, supercooling, and
pressure-melt diagnostics.

## Practical guidance

- Use `absuc` for additive error floors; `reluc` for proportional uncertainty.
- Keep `correlation_length` meaningfully smaller than domain size (as a rule of
  thumb, at least around 3x smaller than domain extent).
- Set `Random.seed!` before building samplers for reproducible runs.
- Start with small `n` while validating model setup; scale up after checks.

## Examples

- `examples/wwfr-simple.jl`: quick start
- `examples/randomly/source-uncertainty-sweep.jl`: source-vs-DEM uncertainty
  comparison and sensitivity sweep

From the `examples/` environment:

```julia
include("wwfr-simple.jl")
include("randomly/source-uncertainty-sweep.jl")
```

See also: [Examples](@ref ExamplesPage), [Subglacially](@ref SubglaciallyGuide), and
[API Reference](@ref APIReference).
