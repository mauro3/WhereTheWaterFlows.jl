using Serialization

@testset "Explicit GRF RNGs" begin
    nx, ny, len = 24, 22, 2.0
    for T in (Float32, Float64), kernel in (WWFR.gaussian_kernel, WWFR.exponential_kernel)
        plain = WWFR.make_grf_sampler(nx, ny, kernel, len; T)
        ker = WWFR.make_kernel(kernel, nx, ny, T(len))
        expected = WWFR.eq18(Xoshiro(73), nx, ny, T(len), ker)
        @test plain(Xoshiro(73)) == expected
        @test plain(Xoshiro(73)) == plain(Xoshiro(73))
        @test plain(Xoshiro(73)) != plain(Xoshiro(74))
        @test eltype(expected) == T
        @test plain(MersenneTwister(73)) == plain(MersenneTwister(73))
        for fftw_plan in (:estimate, :measure)
            Random.seed!(198)
            next_global = rand()
            Random.seed!(198)
            planned = WWFR.make_grf_sampler(nx, ny, kernel, len; T, fftw_plan)
            actual = planned(Xoshiro(73))
            @test rand() == next_global # neither planning nor explicit sampling consumes it
            @test actual == planned(Xoshiro(73))
            @test actual ≈ expected
            @test eltype(actual) == T
        end
        Random.seed!(832)
        legacy = plain()
        Random.seed!(832)
        @test legacy == plain(Random.default_rng())
    end

    field = reshape(collect(1.0:528.0), nx, ny)
    for bounds in ((-Inf, Inf), (2.0, 3.0))
        uc = WWFR.Uncertainty(;abs_bounds=bounds)
        sampler = WWFR.make_sampler(1.0, field, uc)
        rng = Xoshiro(52)
        next_draw = rand(copy(rng))
        @test WWFR.make_field_realization(rng, field, sampler, uc) ==
              field .+ clamp(0.0, bounds...)
        @test rand(rng) == next_draw
    end
    uc = WWFR.Uncertainty(absuc=0.3, reluc=0.2, correlation_length=2.0,
                          abs_bounds=(-0.5, 0.5))
    sampler = WWFR.make_sampler(1.0, field, uc)
    noise = sampler(Xoshiro(73))
    @test WWFR.make_field_realization(Xoshiro(73), field, sampler, uc) ==
          field .+ clamp.(noise .* (0.3 .+ 0.2 .* field), -0.5, 0.5)
end

@testset "Per-field seed table" begin
    names = (:surface, :bed, :pressure)
    seeds = WWFR.make_sample_seeds(1000, names; seed=42)
    @test length(seeds) == 1000
    @test keys(first(seeds)) == names
    @test all(s -> all(v -> v isa UInt64, s), seeds)
    @test seeds == WWFR.make_sample_seeds(1000, names; seed=42)
    @test seeds[1:73] == WWFR.make_sample_seeds(73, names; seed=42)
    @test seeds != WWFR.make_sample_seeds(1000, names; seed=43)
    more = WWFR.make_sample_seeds(1000, (:melt, :pressure, :surface, :bed); seed=42)
    for name in names
        @test getproperty.(seeds, name) == getproperty.(more, name)
    end
    @test seeds[73].surface != seeds[73].bed
    Random.seed!(148)
    next_global = rand()
    Random.seed!(148)
    WWFR.make_sample_seeds(3, names; seed=42)
    @test rand() == next_global
    @test isempty(WWFR.make_sample_seeds(0, names; seed=42))
end

@testset "Seeded routing samples and accumulated provenance" begin
    nx, ny, dx = 24, 22, 100.0
    bed = [1000.0-5i-j for i in 1:nx, j in 1:ny]
    surface = bed .+ 100.0
    ff = fill(0.97, nx, ny)
    source = fill(1e-8, nx, ny)
    baselines = copy.((surface, bed, ff, source))
    groups = [WWF.waterflows(bed).sinks]
    uncertain = WWFR.Uncertainty(absuc=0.1, correlation_length=200.0)
    ff_uc = WWFR.Uncertainty(absuc=0.001, correlation_length=200.0)
    source_uc = WWFR.Uncertainty(reluc=0.1, correlation_length=200.0)
    fixed = WWFR.Uncertainty()
    for kind in (:subaerial, :subglacial)
        build(src_uc; gamma=0.0) = kind == :subaerial ?
            WWFR.make_fns_subaerial(dx, bed, uncertain, source, src_uc, groups) :
            WWFR.make_fns_subglacial(dx, surface, uncertain, bed, uncertain,
                ff, ff_uc, source, src_uc, groups; gamma)
        fns = build(source_uc)
        fields = kind == :subaerial ? (:dem, :source) : (:surface, :bed, :floatfrac, :source)
        seeds = WWFR.make_sample_seeds(1000, fields; seed=42)
        sample73 = fns.sample(;seeds=seeds[73])
        @test sample73 isa WWFR.SeededSample
        @test sample73.seeds == seeds[73]
        @test sample73.rng_type === Xoshiro
        # Direct access to draw1000 needs only its fields, not previous large arrays.
        fns.sample(;seeds=seeds[1000])
        randn(100)
        @test sample73.fields == fns.sample(;seeds=seeds[73]).fields
        @test sample73.fields != fns.sample(;seeds=seeds[1000]).fields
        # Task scheduling and disabled source uncertainty do not shift other fields.
        tasks = [Threads.@spawn fns.sample(;seeds=s) for s in seeds[1:3]]
        @test [fetch(t).fields for t in tasks] == [fns.sample(;seeds=s).fields for s in seeds[1:3]]
        fixed_source = build(fixed).sample(;seeds=seeds[73])
        @test fixed_source.fields[1:end-1] == sample73.fields[1:end-1]
        @test fixed_source.fields[end] == source
        @test build(source_uc;gamma=-0.31).sample(;seeds=seeds[73]).fields == sample73.fields

        ids = [73, 1000, 12]
        selected = seeds[ids]
        aggr = WWFR.map_mc(fns.model, fns.sample, fns.initialize,
            fns.aggregate!, fns.finalize!, length(ids); seeds=selected, sample_ids=ids,
            progressmeter=false)
        @test aggr.sample_ids == ids
        @test aggr.sample_seeds == selected
        @test aggr.rng_type === Xoshiro
        @test aggr.n_samples[] == length(aggr.sample_seeds) == length(ids)
        # Compare each recorded flux and the final maps to independent direct replay.
        manual = fns.initialize()
        for s in aggr.sample_seeds
            fns.aggregate!(manual, fns.model(fns.sample(;seeds=s).fields...))
        end
        fns.finalize!(manual)
        for key in keys(manual)
            a, b = getproperty(aggr, key), getproperty(manual, key)
            if key == :n_samples
                @test a[] == b[]
            elseif key == :catchments
                @test Array(a) == Array(b)
            else
                @test isequal(a, b)
            end
        end
        reversed = WWFR.map_mc(fns.model, fns.sample, fns.initialize,
            fns.aggregate!, fns.finalize!, length(ids); seeds=reverse(selected),
            sample_ids=reverse(ids), progressmeter=false)
        fluxes(a) = kind == :subaerial ? a.catchment_fluxes : a.catchment_fluxes.total
        @test fluxes(reversed) == reverse.(fluxes(aggr))
        @test reversed.sample_ids == reverse(ids)
        buffer = IOBuffer()
        serialize(buffer, aggr)
        seekstart(buffer)
        restored = deserialize(buffer)
        @test restored.sample_ids == ids
        @test restored.sample_seeds == selected
        @test restored.rng_type === Xoshiro
        @test fns.sample(;seeds=restored.sample_seeds[1]).fields == sample73.fields
        legacy = WWFR.map_mc(fns.model, fns.sample, fns.initialize,
            fns.aggregate!, fns.finalize!, 1; progressmeter=false)
        @test !hasproperty(legacy, :sample_seeds)
    end
    @test (surface, bed, ff, source) == baselines
end

@testset "Custom seeded callbacks and failures" begin
    seeds = WWFR.make_sample_seeds(3, (:pressure,); seed=42)
    calls = Ref(0)
    function sample(;seeds)
        calls[] += 1
        WWFR.SeededSample((randn(Xoshiro(seeds.pressure)),), seeds)
    end
    initialize() = (;values=Float64[])
    latest = Ref{Any}(nothing)
    function aggregate!(aggr, value)
        latest[] = aggr
        push!(aggr.values, value)
    end
    aggr = WWFR.map_mc(identity, sample, initialize, aggregate!, identity, 1;
        seeds=seeds[3:3], sample_ids=[1000], progressmeter=false)
    @test calls[] == 1
    @test aggr.sample_ids == [1000]
    @test aggr.sample_seeds == seeds[3:3]
    @test aggr.values == [randn(Xoshiro(seeds[3].pressure))]

    count = Ref(0)
    model(x) = ((count[] += 1) == 3 ? error("model failure") : x)
    @test_throws ErrorException WWFR.map_mc(model, sample, initialize, aggregate!, identity, 3;
        seeds, sample_ids=[71, 72, 73], progressmeter=false)
    @test latest[].sample_ids == [71, 72]
    @test latest[].sample_seeds == seeds[1:2]
    @test length(latest[].values) == 2

    count[] = 0
    function failing_aggregate!(aggr, value)
        latest[] = aggr
        (count[] += 1) == 3 && error("aggregation failure")
        push!(aggr.values, value)
    end
    @test_throws ErrorException WWFR.map_mc(identity, sample, initialize,
        failing_aggregate!, identity, 3; seeds, progressmeter=false)
    @test latest[].sample_ids == [1, 2]
    @test latest[].sample_seeds == seeds[1:2]
    @test length(latest[].values) == 2

    @test_throws DimensionMismatch WWFR.map_mc(identity, sample, initialize,
        aggregate!, identity, 2; seeds, progressmeter=false)
    @test_throws DimensionMismatch WWFR.map_mc(identity, sample, initialize,
        aggregate!, identity, 3; seeds, sample_ids=[1], progressmeter=false)
    @test_throws ArgumentError WWFR.map_mc(identity, sample, initialize,
        aggregate!, identity, 0; progressmeter=false)
    @test_throws TypeError WWFR.map_mc(identity, (;seeds) -> (1.0,), initialize,
        aggregate!, identity, 3; seeds, progressmeter=false)
    @test_throws TypeError WWFR.map_mc(identity, sample, () -> Float64[],
        aggregate!, identity, 3; seeds, progressmeter=false)

    # Do not claim a single RNG/schema if a custom sampler changes it partway through.
    count[] = 0
    function changing_rng(;seeds)
        rng_type = (count[] += 1) == 1 ? Xoshiro : MersenneTwister
        WWFR.SeededSample((randn(rng_type(seeds.pressure)),), seeds; rng_type)
    end
    @test_throws ArgumentError WWFR.map_mc(identity, changing_rng, initialize,
        aggregate!, identity, 3; seeds, progressmeter=false)
    @test length(latest[].values) == length(latest[].sample_seeds) == 1
    count[] = 0
    function changing_schema(;seeds)
        record = (count[] += 1) == 1 ? seeds : (;other=seeds.pressure)
        WWFR.SeededSample((1.0,), record)
    end
    @test_throws ArgumentError WWFR.map_mc(identity, changing_schema, initialize,
        aggregate!, identity, 3; seeds, progressmeter=false)
    @test length(latest[].values) == length(latest[].sample_seeds) == 1
end
