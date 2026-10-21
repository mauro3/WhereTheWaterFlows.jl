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

    area2 = WWF.waterflows(dem, bnd_as_sink=true)[1]
    lakes = demf.>dem
    @test area[.!lakes]==area2[.!lakes]
end

@testset "catchment" begin
    for demfn in [peaks, peaks2, peaks2_nan, peaks2_nan_edge]
        xs, dem = demfn()
        ys = xs

        @test size(dem)==(length(xs), length(ys))
        area, slen, dir, nout, nin, sinks, pits, c, bnds = WWF.waterflows(dem; drain_pits=false);
        @test !isempty(pits)

        for cc in eachindex(pits)
            ij = pits[cc]
            @test catchment(dir, ij) == (c.==cc)
        end

        for (cc,dd) in zip(eachindex(pits), reverse(eachindex(pits)))
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

        @test length(unique(prune_catchments(c, 10))) < length(unique(c))
    end
end
@testset "In-place catchment" begin
    dir = fill(WWF.ind2dir(CartesianIndex(1,0)), 5, 4)
    dir[end,:] .= WWF.SINK
    original = copy(dir)
    first_column = falses(5,4); first_column[:,1] .= true
    last_column = falses(5,4); last_column[:,4] .= true
    first_two = falses(5,4); first_two[:,1:2] .= true
    for mask in (trues(5,4), fill(true,5,4), view(trues(7,6), 2:6, 2:5))
        @test catchment!(mask, dir, CartesianIndex(5,1)) === mask
        @test mask == first_column
        @test catchment!(mask, dir, (5,4)) === mask
        @test mask == last_column # previous catchment was cleared
        seeds = [CartesianIndex(5,1), CartesianIndex(3,1),
                 CartesianIndex(5,1), CartesianIndex(5,4)]
        @test catchment!(mask, dir, seeds) === mask
        @test mask == (first_column .| last_column)
        @test catchment!(mask, dir, CartesianIndices((5:5,1:2))) === mask
        @test mask == first_two
        @test catchment!(mask, dir, CartesianIndex{2}[]) === mask
        @test !any(mask)
    end
    @test catchment(dir, CartesianIndex(5,1)) == first_column
    @test catchment(dir, (5,4)) == last_column
    @test dir == original
    @test_throws DimensionMismatch catchment!(falses(4,4), dir, (5,1))
end

@testset "prune_catchments assigns zero to small catchments" begin
    c = [1 1 2; 1 1 3; 4 4 3]
    original = copy(c)
    @test prune_catchments(c, 2) == [1 1 0; 1 1 2; 3 3 2]
    @test prune_catchments(c, 3) == [1 1 0; 1 1 0; 0 0 0]
    @test prune_catchments(c, 5) == zeros(Int, size(c))
    @test prune_catchments(c, 4) == [1 1 0; 1 1 0; 0 0 0]
    @test prune_catchments(c, 1) == c
    @test prune_catchments([0 -1 2; 4 4 2], 2) == [0 -1 1; 2 2 1]
    @test c == original
    @test prune_catchments([0 2 2; 7 7 7], 0) == [0 1 1; 2 2 2]
    @test prune_catchments([0 2 2; 7 7 7], 3) == [0 0 0; 1 1 1]
    @test prune_catchments(zeros(Int, 2, 3), 2) == zeros(Int, 2, 3)
    @test prune_catchments(fill(-1, 2, 3), 2) == fill(-1, 2, 3)
    @test size(prune_catchments(zeros(Int, 0, 0), 2)) == (0, 0)
end
