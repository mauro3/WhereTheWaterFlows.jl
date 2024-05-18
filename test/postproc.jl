@testset "fill_dem" begin
    # TODO: need more fill tests
    dx = 0.1
    xs = -1.5:dx:1
    ys = -0.5:dx:3.0
    dem = dem1.(xs, ys', withpit=true)
    area, slen, dir, nout, nin, sinks, pits, c, bnds = WWF.waterflows(dem, bnd_as_sink=true)
    demf = WWF.fill_dem(dem, sinks, dir)
    @test sum(demf.-dem) ≈ 2.1499674517313414
    @test sum(demf.-dem .> 0) == 40
    @test all([c[pits[cc]]==cc  for cc=axes(pits)[1]]) # pit in catchment of same color
    @test all([c[sinks[cc]]==cc  for cc=axes(sinks)[1]]) # sink in catchment of same color
end

@testset "catchment" begin
    for demfn in [peaks, peaks2, peaks2_nan, peaks2_nan_edge]
        xs, dem = demfn()
        ys = xs

        @test size(dem)==(length(xs), length(ys))
        area, slen, dir, nout, nin, sinks, pits, c, bnds = WWF.waterflows(dem);

        for cc=1:length(pits)
            ij = pits[cc]
            @test catchment(dir, ij) == (c.==cc)
        end

        for (cc,dd) in zip(1:length(pits), length(pits):-1:1)
            ii, jj = pits[cc], pits[dd]
            @test catchment(dir, [ii,jj]) == ((c.==cc) .| (c.==dd))
        end

        ci = CartesianIndices((2:4,7:9))
        @test catchment(dir, ci) == catchment(dir, collect(ci)[:])

        @test catchment_flux(ones(size(dir)), catchment(dir, ci)) > 0
        @test catchment_flux(ones(size(dir)), c, 5) > 0
        @test catchment_flux(zeros(size(dir)), catchment(dir, ci)) == 0

        sinks = [CartesianIndices((2:4,7:9)), CartesianIndices((5:8,1:3))]
        ss = [collect(sinks[1])[:]; collect(sinks[2])[:]]
        @test (catchments(dir, sinks) .> 0) == catchment(dir, ss)

        @test length(unique(prune_catchments(c, 10; val=0))) < length(unique(c))
    end
end
