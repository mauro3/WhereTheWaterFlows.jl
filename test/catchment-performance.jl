using Test, WhereTheWaterFlows

@testset "Overlapping catchment reductions" begin
    WWF = WhereTheWaterFlows
    WWFS, WWFR = WWF.Subglacially, WWF.Randomly
    values = reshape(1.:153., 17, 9)
    sparse = BitMatrix(values .% 7 .== 0)
    for mask in (falses(17,9), trues(17,9), sparse, Matrix(sparse))
        @test WWF.catchment_flux(values, mask) == sum(values[mask])
    end
    @test WWF.catchment_flux(values, Int.(sparse), 1) == sum(values[sparse])
    bed = [1000. - 10i + j for i in 1:7, j in 1:7]
    surf, source, ff = bed .+ 100, fill(1e-8,7,7), fill(0.97,7,7)
    base = WWFS.waterflows_subglacial(surf, bed, 200., ff, source; gamma=0.)
    sinks = base.routing.sinks
    groups = [sinks, [CartesianIndex(4,4)], vcat(sinks, sinks), CartesianIndex{2}[]]
    for gamma in (0., -0.31)
        out = WWFS.waterflows_subglacial(surf, bed, 200., ff, source; gamma, ctch_sinks=groups)
        reference = WWFS.waterflows_subglacial(surf, bed, 200., ff, source; gamma)
        @test out.routing.c == reference.routing.c
        @test typeof(out.routing.c) == typeof(reference.routing.c)
        a = out.routing.area
        for (i, c) in enumerate(out.sink_catchments.masks)
            @test out.sink_catchments.fluxes.total[i] ≈
                sum(((source .+ a.dissipation_melt_rate .+ a.pressure_melt_rate)*200^2)[c])
            @test out.sink_catchments.fluxes.dissipation[i] ≈ sum((a.dissipation_melt_rate*200^2)[c])
            @test out.sink_catchments.fluxes.pressmelt[i] ≈ sum((a.pressure_melt_rate*200^2)[c])
            @test WWF.catchment_flux(source, c) ≈ sum(source[c])
            layer = zeros(Float16, size(c))
            for _ in 1:1000
                WWFR._accumulate_catchment!(layer, c)
            end
            @test layer == Float16.(c) .* 1000
            layer ./= 1000
            @test layer == c
        end
        @test out.sink_catchments.masks[1] == out.sink_catchments.masks[3]
        @test all(out.sink_catchments.masks[1][out.sink_catchments.masks[2]])
    end
    # Exercise recursive traversal at a modest depth; call-stack size is finite.
    dir = fill(WWF.ind2dir(CartesianIndex(1,0)), 1000, 1)
    dir[end] = WWF.SINK
    @test all(WWF.catchment(dir, CartesianIndex(1000,1)))
    fixed = WWFR.Uncertainty()
    model, sample, reduce! = WWFR.make_fns_subglacial(200., surf, fixed, bed, fixed,
        ff, fixed, source, fixed, groups; gamma=0.)
    out = model(sample()...)
    aggr = WWFR.map_mc(model, sample, reduce!, 2; progressmeter=false)
    @test eltype(aggr.catchments) == Float16
    for i in eachindex(groups)
        @test aggr.catchments[:,:,i] == out[2].sink_catchments.masks[i]
    end
end

@testset "Fused stochastic melt aggregation" begin
    WWFR = WhereTheWaterFlows.Randomly
    n = 64
    bed = [1000. - 10i + j for i in 1:n, j in 1:n]
    surface = [2000. - 5i + 0.2j for i in 1:n, j in 1:n]
    source, floatfrac = fill(1e-8,n,n), fill(0.97,n,n)
    fixed = WWFR.Uncertainty()
    groups = Vector{CartesianIndex{2}}[]
    model, _, reduce! = WWFR.make_fns_subglacial(200., surface, fixed, bed, fixed,
        floatfrac, fixed, source, fixed, groups; gamma=-0.31)
    res = model(surface,bed,floatfrac,source)
    area = res[2].routing.area
    # Exercise cancellation, freezing and inactive-cell NaNs in this reducer fixture.
    area.dissipation_melt_rate[1,1:3] .= [1., -2., NaN]
    area.pressure_melt_rate[1,1:3] .= [-1., 0.5, NaN]
    aggr = reduce!()
    expected = zeros(Float32,n,n)
    for _ in 1:3
        # Preserve the original unfused expression as an independent reference.
        expected .+= area.dissipation_melt_rate + area.pressure_melt_rate
        reduce!(aggr,res)
        @test isequal(aggr.melt_rate,expected)
    end
    # Allow the small lake masks/statistics, not a Float64 melt-rate temporary.
    bytes = @allocated reduce!(aggr,res)
    @test bytes < 4 * length(source)
end
