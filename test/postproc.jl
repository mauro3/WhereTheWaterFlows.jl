@testset "fill_dem" begin
    # TODO: need more fill tests
    dx = 0.1
    xs = -1.5:dx:1
    ys = -0.5:dx:3.0
    dem = dem1.(xs, ys', withpit=true)
    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, bnd_as_pits=false)
    demf = WWF.fill_dem(dem, pits, dir)
    @test sum(demf.-dem) ≈ 2.1499674517313414
    @test sum(demf.-dem .> 0) == 5
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color
end

@testset "catchment" begin
    for demfn in [peaks, peaks2, peaks2_nan, peaks2_nan_edge]
        xs, dem = demfn()
        ys = xs

        @test size(dem)==(length(xs), length(ys))
        area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem);

        for cc=1:length(pits)
            ij = pits[cc]
            @test catchment(dir, ij) == (c.==cc)
        end

        for (cc,dd) in zip(1:length(pits), length(pits):-1:1)
            ii, jj = pits[cc], pits[dd]
            @test catchment(dir, [ii,jj]) == ((c.==cc) .| (c.==dd))
        end
    end
end
