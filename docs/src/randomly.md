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

so kernels operate in grid-cell units.

### Kernel choice

- `gaussian_kernel`: smoother, more diffuse perturbations
- `exponential_kernel`: rougher, shorter-scale structure for same nominal length
- custom kernel: pass any `kernel_fn(nx, ny, len)` compatible function

## What can be configured

### High-level Monte Carlo controls

- `map_mc(model, sample, reduce!, n; progressmeter=true)`
- `n`: number of realizations (`n <= 2047` enforced)
- `progressmeter`: set `false` for quiet batch/CI runs

### Subaerial wrapper (`make_fns_subaerial`)

- uncertain fields: `dem`, `source`
- routing controls: `drain_pits`, `bnd_as_sink`, `nan_as_sink`
- sink groups for aggregation: `ctch_sinks`

### Subglacial wrapper (`make_fns_subglacial`)

- uncertain fields: `surfdem`, `beddem`, `floatfrac`, `source`
- physical controls: `gamma`, `rhow`, `rhoi`
- processing controls: `mask`, `ctch_sinks`, `min_lake_depth`

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
