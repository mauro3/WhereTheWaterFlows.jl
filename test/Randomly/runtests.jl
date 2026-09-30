module Randomly
using Test
using Random
using WhereTheWaterFlows

const WWFR = WhereTheWaterFlows.Randomly
const WWFS = WhereTheWaterFlows.Subglacially
const WWF = WhereTheWaterFlows

const _offset = [0.2, 0.5]

"DEM with a few more maxs and mins. NaN masked edge."
function peaks2_nan_edge(n=100, randfac=0.0)
    coords = range(-pi+_offset[1], pi-_offset[2], length=n)
    dem = sin.(coords) .* cos.(coords') .-
        0.7*(sin.(coords .+ 1) .* cos.(coords')).^8 .+
        0.01 * (coords .+ coords') .+
        randfac * randn(n, n)
    outflow = dem[n÷2, 1]
    dem[1, :] .= NaN
    dem[end, :] .= NaN
    dem[:, 1] .= NaN
    dem[:, end] .= NaN
    dem[n÷2, 1] = outflow
    return coords, dem
end

"Smooth ice-surface DEM"
function icesurf_nan_edge(n=100, randfac=0.0)
    coords = range(-pi+_offset[1], pi-_offset[2], length=n)
    dem = 0.05 * coords .+ (0.05 * coords') .+ 4 .+ randfac * randn(n, n)

    outflow = dem[n÷2, 1]
    dem[1, :] .= NaN
    dem[end, :] .= NaN
    dem[:, 1] .= NaN
    dem[:, end] .= NaN
    dem[n÷2, 1] = outflow
    return coords, dem
end

@testset "GRF" begin
    @testset "find_next_2357" begin
        @test WWFR.find_next_2357(11) == 12
        @test WWFR.find_next_2357(23) == 24
        @test WWFR.find_next_2357(100) == 100
    end

    @testset "kernels and sampler" begin
        nx, ny, len = 61, 63, 5.0
        g_kernel = WWFR.gaussian_kernel(nx, ny, len)
        e_kernel = WWFR.exponential_kernel(nx, ny, len)
        @test maximum(g_kernel) == 1
        @test maximum(e_kernel) == 1

        sa_gaussian = WWFR.make_grf_sampler(128, 128, WWFR.gaussian_kernel, 10.0)
        sa_exponential = WWFR.make_grf_sampler(128, 128, WWFR.exponential_kernel, 10.0)
        @test size(sa_gaussian()) == (128, 128)
        @test size(sa_exponential()) == (128, 128)
    end
end

@testset "Zero-uncertainty sampling" begin
    field = Float32[1 2; NaN -3]
    original = copy(field)
    unused_kernel(args...) = error("Zero uncertainty must not build a GRF kernel")
    for (absuc, reluc) in ((0, 0), (zeros(Float32, 2, 2), 0),
                           (0, zeros(Float32, 2, 2)),
                           (zeros(Float32, 2, 2), zeros(Float32, 2, 2)))
        for bounds in ((-Inf, Inf), (2.0, 3.0), (-3.0, -2.0))
            uc = WWFR.Uncertainty(; absuc, reluc, abs_bounds=bounds,
                                   covariance_fn=unused_kernel)
            sampler = WWFR.make_sampler(1.0, field, uc)
            @test iszero(sampler())
            result = WWFR.make_field_realization(field, sampler, uc)
            @test isequal(result, field .+ Float32(clamp(0.0, bounds...)))
            @test eltype(result) == Float32
            @test result !== field
            result[1, 1] = 100
            @test isequal(field, original)
            @test isequal(WWFR.make_field_realization(field, sampler, uc),
                          field .+ Float32(clamp(0.0, bounds...)))
        end
    end

    # Neither constructing nor invoking a zero-uncertainty sampler consumes RNG draws.
    Random.seed!(42)
    expected_draw = rand()
    Random.seed!(42)
    uc = WWFR.Uncertainty()
    sampler = WWFR.make_sampler(1.0, field, uc)
    WWFR.make_field_realization(field, sampler, uc)
    @test rand() == expected_draw

    # A nonzero entry in either uncertainty field still takes the GRF path.
    nonzero = Float32[0 0; 0 1]
    for (absuc, reluc) in ((nonzero, 0), (0, nonzero))
        uc = WWFR.Uncertainty(; absuc, reluc, correlation_length=0.1,
                               covariance_fn=unused_kernel)
        @test_throws ErrorException WWFR.make_sampler(1.0, field, uc)
    end
end

@testset "Stochastic Subglacial" begin
    n = 100
    x, beddem = peaks2_nan_edge(n)
    x, surfdem = icesurf_nan_edge(n)
    x = x .* 10000
    dx = step(x)
    beddem .*= 20
    surfdem .*= 100
    y = x
    mask = fill!(similar(surfdem, Bool), true)

    source = fill!(similar(surfdem), 1) / dx^2
    floatfrac = fill!(similar(surfdem), 1)
    ctch_sinks = [CartesianIndices((1:10, 1:length(y)))[:], CartesianIndices((1:10, 1:length(y)÷2))[:]]

    surfdem_uc = WWFR.Uncertainty(absuc=0, reluc=0.0, correlation_length=1e4)
    beddem_uc = WWFR.Uncertainty(absuc=0, reluc=0.0, correlation_length=1e3)
    floatfrac_uc = WWFR.Uncertainty(absuc=0, reluc=0.0, correlation_length=1e3)
    source_uc = WWFR.Uncertainty()

    model, sample, reduce! = WWFR.make_fns_subglacial(dx,
                                                      surfdem, surfdem_uc,
                                                      beddem, beddem_uc,
                                                      floatfrac, floatfrac_uc,
                                                      source, source_uc,
                                                      ctch_sinks;
                                                      mask)

    # reluc==0 should be deterministic and match direct WWFS call
    input, output = model(sample()...)
    out_ref = WWFS.waterflows_subglacial(surfdem, beddem, dx, floatfrac, source, mask;
                                         ctch_sinks=ctch_sinks,
                                         gamma=0,
                                         rhow=WWFS.RHOW,
                                         rhoi=WWFS.RHOI,
                                         drain_pits=true,
                                         bnd_as_sink=true,
                                         nan_as_sink=true)
    @test output.routing.dir == out_ref.routing.dir
    @test output.routing.c == out_ref.routing.c
    @test output.pressmelt.kappas == out_ref.pressmelt.kappas

    aggr = WWFR.map_mc(model, sample, reduce!, 3)
    @test size(aggr.areas_total) == size(surfdem)
    @test size(aggr.lake_depth_free_surface) == size(surfdem)
    @test length(aggr.catchment_fluxes.total) == length(ctch_sinks)
    @test aggr.n_samples[] == 3
end

@testset "Subaerial catchment flux mass balance" begin
    dem = [Float64(i + j) for i in 1:5, j in 1:5]
    dx = 3.0
    baseline = ones(size(dem))
    sinks = WWF.waterflows(dem).sinks
    # Including an upstream point in the all-outlets group must not count it twice.
    groups = [sinks, [sinks; CartesianIndex(3, 3)], [CartesianIndex(3, 3)]]
    model, _, reduce! = WWFR.make_fns_subaerial(
        dx, dem, WWFR.Uncertainty(), baseline, WWFR.Uncertainty(), groups)
    aggr = reduce!()
    for scale in (1.0, 2.0)
        source = scale .* reshape(collect(1.0:25.0), size(dem))
        res = model(dem, source)
        reduce!(aggr, res)
        out = res[2]
        total = sum(source) * dx^2
        @test aggr.catchment_fluxes[1][end] ≈ total
        @test aggr.catchment_fluxes[2][end] ≈ total
        @test aggr.catchment_fluxes[1][end] ≈ sum(out.area[sinks])
        @test aggr.catchment_fluxes[3][end] ≈ out.area[3, 3]
    end
    @test aggr.catchment_fluxes[1][2] ≈ 2 * aggr.catchment_fluxes[1][1]
end

@testset "Stochastic Subaerial" begin
    n = 60
    _, dem = peaks2_nan_edge(n)
    dx = 100.0
    source = fill!(similar(dem), 1.0) / dx^2
    ctch_sinks = [CartesianIndices((2:2, 2:n-1))[:]]

    dem_uc = WWFR.Uncertainty(absuc=0, reluc=0.0, correlation_length=10.0)
    source_uc = WWFR.Uncertainty()

    model, sample, reduce! = WWFR.make_fns_subaerial(dx,
                                                     dem, dem_uc,
                                                     source, source_uc,
                                                     ctch_sinks)

    input, output = model(sample()...)
    out_ref = WWF.waterflows(dem, source .* dx^2; drain_pits=true, bnd_as_sink=true, nan_as_sink=true)

    @test output.dir == out_ref.dir
    @test output.c == out_ref.c
    @test output.slen == out_ref.slen

    aggr = WWFR.map_mc(model, sample, reduce!, 3)
    @test size(aggr.areas_total) == size(dem)
    @test size(aggr.catchments, 3) == length(ctch_sinks)
    @test aggr.n_samples[] == 3
end
end
