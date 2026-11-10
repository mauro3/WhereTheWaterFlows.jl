module Randomly

using ..WhereTheWaterFlows
using ..Subglacially

using ProgressMeter

include("grf.jl")
include("blocksparse.jl")

export Uncertainty, map_mc, make_fns_subaerial, make_fns_subglacial,
       gaussian_kernel, exponential_kernel

"""
Type to hold the uncertainty information for a field.

- absuc=0 -- absolute uncertainty as std (array or scalar)
- reluc=0 -- relative (to local field value) uncertainty as std (array or scalar)
- correlation_length=1 -- ditto (scalar)
- abs_bounds=(-Inf,Inf) -- values of the GRF are constrained to lie within these bounds (simple crop,
                            after being scaled with absuc and reluc)

When both uncertainties are zero (including all-zero arrays), sampling skips
GRF setup and FFTs and consumes no random draws. Perturbation bounds still apply.
"""
struct Uncertainty{A,R,C}
    absuc::A
    reluc::R
    correlation_length::Float64
    covariance_fn::C
    abs_bounds::Tuple{Float64,Float64}
end

function Uncertainty(; absuc=0,
                       reluc=0,
                       correlation_length=1,
                       covariance_fn=gaussian_kernel,
                       abs_bounds=(-Inf, Inf))
    return Uncertainty(absuc, reluc, Float64(correlation_length),
                       covariance_fn, Float64.(abs_bounds))
end

function make_sampler(dx, field, uncert)
    if iszero(uncert.absuc) && iszero(uncert.reluc)
        # Keep the usual scaling/bounds path, but avoid kernel setup and FFTs.
        return () -> fill!(similar(field), zero(eltype(field)))
    end
    (;correlation_length, covariance_fn) = uncert
    len = correlation_length/dx
    nx, ny = size(field)
    T = eltype(field)

    return make_grf_sampler(nx, ny, covariance_fn, len; T, fftw_plan=nothing)
end

"""
    make_field_realization(field, sampler, uncert::Uncertainty)

Make a realization of a field.
"""
function make_field_realization(field, sampler, uncert::Uncertainty)
    (;absuc, reluc, abs_bounds) = uncert
    s = sampler()

    # scale amplitude with the error
    @. s = s*absuc + field*s*reluc

    if abs_bounds!=(-Inf, Inf)
        s .= clamp.(s, abs_bounds...)
    end

    # add to field
    return field .+ s
end

"""
    map_mc(model, sample, initialize, aggregate!, finalize!, n; progressmeter=true)

Run Monte Carlo sampling and aggregate model outputs in place.

Args
- model   -- the model function, runs with `output = model(sample()...)`
- sample  -- `input = sample()` returns a new sample to use as input to model
- initialize -- called once to create aggregation storage (no arguments)
- aggregate! -- called as `aggregate!(aggr, model(sample()...))` for each of the `n` realizations
- finalize! -- called as `finalize!(aggr)` once at the end and updates `aggr` in-place.

The sample count must be at most 2048 because catchment frequencies use Float16
counts. Set `progressmeter=false` to suppress progress updates.
"""
function map_mc(model, sample, initialize, aggregate!, finalize!, n; progressmeter=true)
    n>2048 && error("""Cannot take more than 2048 samples as otherwise
                       the aggregation for the catchment probabilities does not work anymore
                       as Float16 are used.""")

    lockobj = Threads.ReentrantLock()
    # initialize aggregate storage
    aggr = initialize()
    p = Progress(n)
    # Threading issues:
    # - excessive memory usage (mem-leak? Maybe fixed now with hand-rolled GRF)
    for _=1:n
        input = sample()
        res = model(input...)
        lock(lockobj) do
            aggregate!(aggr, res)
        end
        progressmeter && next!(p)
    end
    # finalize aggregate
    finalize!(aggr)
    return aggr
end

"""
    make_fns_subglacial(dx,
                          surfdem, surfdem_uc,
                          beddem, beddem_uc,
                          floatfrac, floatfrac_uc,
                          source, source_uc,
                          ctch_sinks;
                          mask::AbstractMatrix=fill!(similar(surfdem, Bool), true),
                          gamma=0.0,
                          min_lake_depth=10.0, # default min lake depth under which value are not aggregated
                          rhow=RHOW, rhoi=RHOI,
                          catchment_storage=:blocksparse,
                          catchment_blocksize=_CATCHMENT_BLOCKSIZE)

Build a named tuple of callbacks for stochastic subglacial routing with `map_mc`.

Args
- `dx` -- grid spacing
- `surfdem, beddem, floatfrac, source` -- baseline fields used for sampling and routing
- `floatfrac` must be an array with the same size as `surfdem`.
  A scalar is not accepted here; use `fill(value, size(surfdem)...)` for a
  spatially uniform flotation field.
- `surfdem_uc, beddem_uc, floatfrac_uc, source_uc` -- corresponding `Uncertainty` objects
- `ctch_sinks` -- sink groups used for catchment masks/flux aggregation

Kwargs
- `mask` -- routing mask passed to `Subglacially.waterflows_subglacial`
- `gamma` -- pressure-melting deflection parameter; defaults to `0.0` when omitted. Otherwise
             set to WWF.GAMMA for a standard value.
- `min_lake_depth` -- strict depth threshold for lake occurrence and volume; volumes
  sum full qualifying depths times `dx^2`, not excess depths above the threshold
- `rhow`, `rhoi` -- water and ice density overrides
- `catchment_storage` -- `:blocksparse` (default) for lazy-zero tiled storage,
  or `:dense` for a conventional Float16 array. Both retain `(x, y, group)` indexing.
  Counts and final probabilities use the same storage; other outputs are unchanged.
- `catchment_blocksize` -- square tile size for `:blocksparse` storage; defaults to
  50. It changes memory use and update performance, not results, and is ignored for
  `:dense` storage.

Return `(; model, sample, initialize, aggregate!, finalize!)`:
- `model(surf, bed, floatfrac, source)` returns
  `(;surf, bed, dx, floatfrac, source), waterflows_subglacial(...)`.
- `sample()` returns `surf, bed, float, source`.
- `initialize()` creates fresh aggregate storage.
- `aggregate!(aggr, out)` accumulates one model output in place.
- `finalize!(aggr)` normalizes accumulated maps, leaving per-sample vectors unchanged.
"""
function make_fns_subglacial(dx,
                              surfdem::AbstractMatrix, surfdem_uc,
                              beddem::AbstractMatrix, beddem_uc,
                              floatfrac::AbstractMatrix, floatfrac_uc,
                              source::AbstractMatrix, source_uc,
                              ctch_sinks;
                              mask::AbstractMatrix=fill!(similar(surfdem, Bool), true),
                              gamma=0.0,
                              min_lake_depth=10.0, # default min lake depth under which value are not aggregated
                              rhow=Subglacially.RHOW, rhoi=Subglacially.RHOI,
                              catchment_storage=:blocksparse,
                              catchment_blocksize=_CATCHMENT_BLOCKSIZE)

    size(floatfrac) == size(surfdem) || throw(DimensionMismatch(
        "floatfrac must have the same size as surfdem"))

    model(surf, bed, floatfrac, source) = ((;surf, bed, dx, floatfrac, source),
                                             Subglacially.waterflows_subglacial(surf, bed, dx, floatfrac, source, mask;
                                                                                gamma, drain_pits=true,
                                                                                bnd_as_sink=true,
                                                                                nan_as_sink=true,
                                                                                rhow, rhoi, ctch_sinks))

    sample = let
        # pre-compute the GRF samplers
        surfdem_grf_sampler = make_sampler(dx, surfdem, surfdem_uc)
        beddem_grf_sampler = make_sampler(dx, beddem, beddem_uc)
        floatfrac_grf_sampler = make_sampler(dx, floatfrac, floatfrac_uc)
        source_grf_sampler = make_sampler(dx, source, source_uc)
        function ()
            surf = make_field_realization(surfdem, surfdem_grf_sampler, surfdem_uc)
            bed = make_field_realization(beddem, beddem_grf_sampler, beddem_uc)
            float = make_field_realization(floatfrac, floatfrac_grf_sampler, floatfrac_uc)
            src = make_field_realization(source, source_grf_sampler, source_uc)
            return surf, bed, float, src
        end
    end

    # initialize storage
    function initialize()
        # TODO maybe use Float16 (but that only goes to 65000!)
        sz = size(surfdem)
        (;areas_total = zeros(Float32,sz),
          areas_extra = zeros(Float32,sz),
          melt_rate = zeros(Float32,sz),
          lake_depth_fixed_surface = zeros(Float32,sz),
          lake_mask_fixed_surface = zeros(Float32,sz),
          lake_depth_free_surface = zeros(Float32,sz),
          lake_mask_free_surface = zeros(Float32,sz),
          sc_locs = zeros(Float32,sz),
          kappas = zeros(Float32,sz),
          catchments = _init_catchments(sz, length(ctch_sinks), catchment_storage;
                                        blocksize=catchment_blocksize),
          catchment_fluxes = (total=[Float32[] for s in ctch_sinks],
                              dissipation=[Float32[] for s in ctch_sinks],
                              pressmelt=[Float32[] for s in ctch_sinks]),
          n_samples = Ref(0),
          lake_vol_fixed_surface = Float32[],
          lake_vol_free_surface = Float32[],
          )
    end

    # aggregate values
    function aggregate!(aggr, res)
        # unpack
        input, output = res
        (;area) = output.routing
        (;sc_locs, kappas) = output.pressmelt
        (;depth_fixed_surface, depth_free_surface) = output.lakes
        (;masks, fluxes) = output.sink_catchments

        aggr.areas_total .+= area.total
        aggr.areas_extra .+= area.extra
        aggr.melt_rate .+= area.dissipation_melt_rate .+ area.pressure_melt_rate
        aggr.lake_depth_fixed_surface .+= depth_fixed_surface
        aggr.lake_mask_fixed_surface .+= depth_fixed_surface.>min_lake_depth

        aggr.lake_depth_free_surface .+= depth_free_surface
        aggr.lake_mask_free_surface .+= depth_free_surface.>min_lake_depth

        push!(aggr.lake_vol_fixed_surface,
              sum(Float32.(depth_fixed_surface[depth_fixed_surface .> min_lake_depth])) * dx^2)
        push!(aggr.lake_vol_free_surface,
              sum(Float32.(depth_free_surface[depth_free_surface .> min_lake_depth])) * dx^2)

        aggr.sc_locs .+= sc_locs
        aggr.kappas .+= kappas

        for i=1:length(fluxes.total)
            push!(aggr.catchment_fluxes.total[i], fluxes.total[i])
            push!(aggr.catchment_fluxes.dissipation[i], fluxes.dissipation[i])
            push!(aggr.catchment_fluxes.pressmelt[i], fluxes.pressmelt[i])
            _accumulate_catchment!(aggr.catchments, masks[i], i)
        end

        aggr.n_samples[] += 1
        return aggr
    end

    # Finalize by taking the mean (except for the catchment_fluxes, which are separate records
    # for each forward model run)
    function finalize!(aggr)
        n = aggr.n_samples[]
        aggr.areas_total ./= n
        aggr.areas_extra ./= n
        aggr.melt_rate ./= n
        aggr.lake_depth_fixed_surface ./= n
        aggr.lake_mask_fixed_surface ./= n

        aggr.lake_depth_free_surface ./= n
        aggr.lake_mask_free_surface ./= n

        aggr.sc_locs ./= n
        aggr.kappas ./= n

        _normalize_catchments!(aggr.catchments, n)

        return aggr
    end

    return (; model, sample, initialize, aggregate!, finalize!)
end

"""
    make_fns_subaerial(dx,
                        dem, dem_uc,
                        source, source_uc,
                        ctch_sinks;
                        drain_pits=true,
                        bnd_as_sink=true,
                        nan_as_sink=true,
                        catchment_storage=:blocksparse,
                        catchment_blocksize=_CATCHMENT_BLOCKSIZE)

Build a named tuple of callbacks for stochastic subaerial routing with `map_mc`.

Args
- `dx` -- grid spacing
- `dem, source` -- baseline elevation and source fields
- `dem_uc, source_uc` -- corresponding `Uncertainty` objects
- `ctch_sinks` -- sink groups used for catchment masks/flux aggregation

Kwargs
- `drain_pits`, `bnd_as_sink`, `nan_as_sink` -- forwarded to `WhereTheWaterFlows.waterflows`
- `catchment_storage` -- `:blocksparse` (default) for lazy-zero tiled storage,
  or `:dense` for a conventional Float16 array. Both retain `(x, y, group)` indexing.
  Counts and final probabilities use the same storage; other outputs are unchanged.
- `catchment_blocksize` -- square tile size for `:blocksparse` storage; defaults to
  50. It changes memory use and update performance, not results, and is ignored for
  `:dense` storage.

Return `(; model, sample, initialize, aggregate!, finalize!)`:
- `model(dem_, source_)` returns `(; dem, dx, source), waterflows(...)`.
- `sample()` returns one sampled `(dem_, source_)`.
- `initialize()` creates fresh aggregate storage.
- `aggregate!(aggr, out)` accumulates one model output in place.
- `finalize!(aggr)` normalizes accumulated maps, leaving per-sample vectors unchanged.
"""
function make_fns_subaerial(dx,
                            dem, dem_uc,
                            source, source_uc,
                            ctch_sinks;
                            drain_pits=true,
                            bnd_as_sink=true,
                            nan_as_sink=true,
                            catchment_storage=:blocksparse,
                            catchment_blocksize=_CATCHMENT_BLOCKSIZE)

    model(dem_, source_) = ((;dem=dem_, dx, source=source_),
                            WhereTheWaterFlows.waterflows(dem_, source_ .* dx^2;
                                                          drain_pits, bnd_as_sink, nan_as_sink))

    sample = let
        dem_sampler = make_sampler(dx, dem, dem_uc)
        source_sampler = make_sampler(dx, source, source_uc)
        function ()
            dem_ = make_field_realization(dem, dem_sampler, dem_uc)
            src_ = make_field_realization(source, source_sampler, source_uc)
            return dem_, src_
        end
    end

    # initialize storage
    function initialize()
        sz = size(dem)
        return (;areas_total = zeros(Float32, sz),
                stream_length = zeros(Float32, sz),
                catchments = _init_catchments(sz, length(ctch_sinks), catchment_storage;
                                              blocksize=catchment_blocksize),
                catchment_fluxes = [Float32[] for _ in ctch_sinks],
                n_samples = Ref(0))
    end

    # aggregate values
    catchment_mask = similar(BitArray, axes(dem)) # scratch storage to avoid allocations
    function aggregate!(aggr, res)
        (input, output) = res
        aggr.areas_total .+= output.area
        aggr.stream_length .+= output.slen

        for i=1:length(ctch_sinks)
            WhereTheWaterFlows.catchment!(catchment_mask, output.dir, ctch_sinks[i])
            _accumulate_catchment!(aggr.catchments, catchment_mask, i)
            flux = WhereTheWaterFlows.catchment_flux(input.source, catchment_mask) * input.dx^2
            push!(aggr.catchment_fluxes[i], flux)
        end

        aggr.n_samples[] += 1
        return aggr
    end

    # finalize output
    function finalize!(aggr)
        n = aggr.n_samples[]
        aggr.areas_total ./= n
        aggr.stream_length ./= n
        _normalize_catchments!(aggr.catchments, n)
        return aggr
    end

    return (; model, sample, initialize, aggregate!, finalize!)
end

end
