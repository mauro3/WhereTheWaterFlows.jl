module Randomly

using ..WhereTheWaterFlows
using ..Subglacially

using ProgressMeter

include("grf.jl")

export Uncertainty, map_mc, make_fns_subaerial, make_fns_subglacial
export find_next_2357, pad, depad, distance_squared, gaussian_kernel, exponential_kernel, make_kernel, make_grf_sampler

"""
Type to hold the uncertainty information for a field.

- absuc=0 -- absolute uncertainty as std (array or scalar)
- reluc=0 -- relative (to local field value) uncertainty as std (array or scalar)
- correlation_length=1 -- ditto (scalar)
- abs_bounds=(-Inf,Inf) -- values of the GRF are constrained to lie within these bounds (simple crop,
                           after being scaled with absuc and reluc)
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
    map_mc(model, sample, reduce!, n)

A general Monte Carlo `map` function

Args
- model   -- the model function, runs with `output = model(sample()...)`
- sample  -- `input = sample()` returns a new sample to use as input to model
- reduce! -- aggregates the model output with `reduce!(aggr, output)`.  Also initializes
                the aggregate storage with `aggr = reduce!()` and finalizes it
                with `reduce!(aggr)`.
- n -- number of samples to take
       (must be <= 2048; Float16 aggregation for catchment frequencies loses
         precision for larger counts)

Thread saftey: only `reduce!` is allowed to not be thread-safe.
"""
function map_mc(model, sample, reduce!, n; progressmeter=true)
    n>2048 && error("""Cannot take more than 2048 samples as otherwise
                       the aggregation for the catchment probabilities does not work anymore
                       as Float16 are used.""")

    lockobj = Threads.ReentrantLock()
    # initialize aggregate! storage
    aggr = reduce!()
    p = Progress(n)
    # Threading issues:
    # - excessive memory usage (mem-leak? Maybe fixed now with hand-rolled GRF)
    for _=1:n
        input = sample()
        res = model(input...)
        lock(lockobj) do
            reduce!(aggr, res)
        end
        progressmeter && next!(p)
    end
    # finalize aggregate
    reduce!(aggr)
    return aggr
end

"""
    make_fns_subglacial(dx,
                          surfdem, surfdem_uc,
                          beddem, beddem_uc,
                          floatfrac, floatfrac_uc,
                          source, source_uc,
                          ctch_sinks;
                          mask=mask::AbstractMatrix=fill!(similar(surfdem, Bool), true),
                          gamma=0.0,
                          min_lake_depth=10.0, # default min lake depth under which value are not aggregated
                          rhow=ROHW, rhoi=RHOW)

Build `model`, `sample`, and `reduce!` functions for stochastic subglacial routing
which are then used in `map_mc`.

Args
- `dx` -- grid spacing
- `surfdem, beddem, floatfrac, source` -- baseline fields used for sampling and routing
- `surfdem_uc, beddem_uc, floatfrac_uc, source_uc` -- corresponding `Uncertainty` objects
- `ctch_sinks` -- sink groups used for catchment masks/flux aggregation

Kwargs
- `mask` -- routing mask passed to `Subglacially.waterflows_subglacial`
- `gamma` -- pressure-melting deflection parameter; defaults to `0.0` when omitted. Otherwise
             set to WWF.GAMMA for a standard value.
- `min_lake_depth` -- threshold used for lake-occurrence masks and lake-volume vectors
- `rhow`, `rhoi` -- water and ice density overrides

Return:
- model(surf, bed, floatfrac, source) -> runs one model realisation returning
  `(;surf, bed, dx, floatfrac, source), waterflows_subglacial(...)`
- sample() -> surf, bed, float, source
- reduce! -> three-method reduction function (`reduce!()`, `reduce!(aggr, out)`, `reduce!(aggr)`).
"""
function make_fns_subglacial(dx,
                              surfdem, surfdem_uc,
                              beddem, beddem_uc,
                              floatfrac, floatfrac_uc,
                              source, source_uc,
                              ctch_sinks;
                              mask=mask::AbstractMatrix=fill!(similar(surfdem, Bool), true),
                              gamma=0.0,
                              min_lake_depth=10.0, # default min lake depth under which value are not aggregated
                              rhow=Subglacially.RHOW, rhoi=Subglacially.RHOI)

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
    function reduce!()
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
          catchments = zeros(Float16,sz...,length(ctch_sinks)),
          catchment_fluxes = (total=[Float32[] for s in ctch_sinks],
                              dissipation=[Float32[] for s in ctch_sinks],
                              pressmelt=[Float32[] for s in ctch_sinks]),
          n_samples = Ref(0),
          lake_vol_fixed_surface = Float32[],
          lake_vol_free_surface = Float32[],
          )
    end

    # aggregate values
    function reduce!(aggr, res)
        # unpack
        input, output = res
        (;area) = output.routing
        (;sc_locs, kappas) = output.pressmelt
        (;depth_fixed_surface, depth_free_surface) = output.lakes
        (;masks, fluxes) = output.sink_catchments

        aggr.areas_total .+= area.total
        aggr.areas_extra .+= area.extra
        aggr.melt_rate .+= area.dissipation_melt_rate + area.pressure_melt_rate
        aggr.lake_depth_fixed_surface .+= depth_fixed_surface
        aggr.lake_mask_fixed_surface .+= depth_fixed_surface.>min_lake_depth

        aggr.lake_depth_free_surface .+= depth_free_surface
        aggr.lake_mask_free_surface .+= depth_free_surface.>min_lake_depth

        push!(aggr.lake_vol_fixed_surface, sum(Float32.(depth_fixed_surface[depth_fixed_surface .> min_lake_depth])))
        push!(aggr.lake_vol_free_surface, sum(Float32.(depth_free_surface[depth_free_surface .> min_lake_depth])))

        aggr.sc_locs .+= sc_locs
        aggr.kappas .+= kappas

        for i=1:length(fluxes.total)
            push!(aggr.catchment_fluxes.total[i], fluxes.total[i])
            push!(aggr.catchment_fluxes.dissipation[i], fluxes.dissipation[i])
            push!(aggr.catchment_fluxes.pressmelt[i], fluxes.pressmelt[i])
            aggr.catchments[:,:,i] .+= masks[i]
        end

        aggr.n_samples[] += 1
        return aggr
    end

    # Finalize by taking the mean (except for the catchment_fluxes, which are separate records
    # for each forward model run)
    function reduce!(aggr)
        n = aggr.n_samples[]
        aggr.areas_total ./= n
        aggr.lake_depth_fixed_surface ./= n
        aggr.lake_mask_fixed_surface ./= n

        aggr.lake_depth_free_surface ./= n
        aggr.lake_mask_free_surface ./= n

        aggr.sc_locs ./= n
        aggr.kappas ./= n

        aggr.catchments ./= n

        return aggr
    end

    return model, sample, reduce!
end

"""
    make_fns_subaerial(dx,
                        dem, dem_uc,
                        source, source_uc,
                        ctch_sinks;
                        drain_pits=true,
                        bnd_as_sink=true,
                        nan_as_sink=true)

Build `model`, `sample`, and `reduce!` functions for stochastic subaerial routing.

Args
- `dx` -- grid spacing
- `dem, source` -- baseline elevation and source fields
- `dem_uc, source_uc` -- corresponding `Uncertainty` objects
- `ctch_sinks` -- sink groups used for catchment masks/flux aggregation

Kwargs
- `drain_pits`, `bnd_as_sink`, `nan_as_sink` -- forwarded to `WhereTheWaterFlows.waterflows`

Return
- model(dem_, source_) -> returns `(; dem, dx, source), waterflows(...)`
- sample() -> returns one sampled `(dem_, source_)`
- reduce! -> three-method reduction function (`reduce!()`, `reduce!(aggr, out)`, `reduce!(aggr)`).
"""
function make_fns_subaerial(dx,
                            dem, dem_uc,
                            source, source_uc,
                            ctch_sinks;
                            drain_pits=true,
                            bnd_as_sink=true,
                            nan_as_sink=true)
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

    function reduce_step(aggr, res; init=false, sz=nothing)
        if init
            return (;areas_total = zeros(Float32, sz),
                    stream_length = zeros(Float32, sz),
                    catchments = zeros(Float16, sz..., length(ctch_sinks)),
                    catchment_fluxes = [Float32[] for _ in ctch_sinks],
                    n_samples = Ref(0))
        end

        (_, output) = res
        aggr.areas_total .+= output.area ./ dx^2
        aggr.stream_length .+= output.slen

        for i=1:length(ctch_sinks)
            c = WhereTheWaterFlows.catchment(output.dir, ctch_sinks[i])
            aggr.catchments[:, :, i] .+= c
            push!(aggr.catchment_fluxes[i], WhereTheWaterFlows.catchment_flux(output.area, c))
        end

        aggr.n_samples[] += 1
        return aggr
    end

    function reduce_finalize(aggr)
        n = aggr.n_samples[]
        aggr.areas_total ./= n
        aggr.stream_length ./= n
        aggr.catchments ./= n
        return aggr
    end

    function reduce!()
        sz = size(dem)
        return reduce_step(nothing, nothing; init=true, sz)
    end
    function reduce!(aggr, res)
        return reduce_step(aggr, res)
    end
    function reduce!(aggr)
        return reduce_finalize(aggr)
    end

    return model, sample, reduce!
end

end
