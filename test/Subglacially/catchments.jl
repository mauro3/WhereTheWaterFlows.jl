@testset "Catchment sink connectivity" begin
    C(i,j) = CartesianIndex(i,j)
    components = WWFS.catchment_sink_components
    repair! = WWFS.connect_catchment_sinks!
    @test isempty(components(CartesianIndex{2}[]))
    @test components([C(1,2),C(2,1)]) == [[C(1,2),C(2,1)]]
    @test length(components([C(1,1),C(2,2),C(3,3)])) == 1
    ring = [C(i,j) for i in 1:3 for j in 1:3 if (i,j)!=(2,2)]
    @test length(components(ring)) == 1
    branch = vcat([C(i,2) for i in 1:5],[C(3,3),C(3,4),C(3,5)])
    @test length(components(branch)) == 1
    @test repair!(Vector{Vector{CartesianIndex{2}}}()) == 0
    @test repair!([CartesianIndex{2}[],CartesianIndex{2}[]]) == 0
    connected = [branch,CartesianIndex{2}[]]
    before = deepcopy(connected)
    @test repair!(connected) == 0
    @test connected == before

    # Repair a small fragment before considering the larger fragment beyond it.
    # The latter now reaches its original group, so only one label needs changing.
    groups = [[C(i,2) for i in (1,2,3,7)], [C(i,2) for i in (4,5,6,8,9)], CartesianIndex{2}[]]
    original = deepcopy(groups)
    @test repair!(groups) == 1
    @test groups == [[C(i,2) for i in 1:3],[C(i,2) for i in 4:9],CartesianIndex{2}[]]
    @test Set(vcat(groups...)) == Set(vcat(original...))
    @test sum(length,groups) == sum(length,original)
    @test all(g -> length(components(g)) <= 1,groups)
    @test isempty(groups[3])
    @test repair!(groups) == 0
    reversed = reverse.(original)
    @test repair!(reversed) == 1
    @test reversed == groups

    # A small fragment may need a later, larger bridge before it becomes reachable.
    chain = [[C(i,2) for i in (1,2,3,4,8,9)], [C(i,2) for i in (5,6,7,12,13,14,15,16,17)]]
    @test repair!(chain) == 3
    @test chain == [[C(i,2) for i in 1:9],[C(i,2) for i in 12:17]]

    # Prefer most contacts; break equal-contact ties by group index.
    tie = [[C(1,1),C(1,2),C(5,5)],[C(3,5),C(4,5)],[C(6,5),C(7,5)]]
    majority = deepcopy(tie)
    push!(majority[3],C(6,6))
    @test repair!(tie) == 1
    @test C(5,5) in tie[2]
    @test repair!(majority) == 1
    @test C(5,5) in majority[3]

    # A fragment can attach to another group's main component diagonally.
    diagonal = [[C(1,1),C(1,2),C(1,3),C(5,5)], [C(6,6),C(6,7)]]
    @test repair!(diagonal) == 1
    @test C(5,5) in diagonal[2]
    @test all(g -> length(components(g)) == 1,diagonal)

    # No fabricated links across genuine gaps; failed repairs leave input unchanged.
    isolated = [[C(1,1),C(1,2),C(20,20)],[C(7,7)]]
    before = deepcopy(isolated)
    @test_throws ArgumentError repair!(isolated)
    @test isolated == before
    partial = deepcopy(original)
    push!(partial[1],C(20,20))
    before = deepcopy(partial)
    @test_throws ArgumentError repair!(partial)
    @test partial == before
    duplicate = [[C(1,1)],[C(1,1)]]
    @test_throws ArgumentError repair!(duplicate)

    # Integration: a polygon cut catches a small nub on the other group's coast.
    sinks = vcat([C(2,j) for j in 2:7],[C(i,7) for i in 3:6],
                 [C(6,j) for j in 8:12],[C(4,10),C(5,10)])
    mask = falses(14,14);mask[sinks] .= true
    polygons = [[(0.,0.),(4.5,0.),(4.5,14.),(0.,14.)],
                [(4.5,0.),(14.,0.),(14.,14.),(4.5,14.)]]
    raw = WWFS.catchment_sinks(polygons,mask,1.:14.,1.:14.,sinks;
        D9_iters=0,D25_iter=0,connect=false)
    @test length(components(raw[1])) == 2
    fixed = WWFS.catchment_sinks(polygons,mask,1.:14.,1.:14.,sinks;
        D9_iters=0,D25_iter=0)
    @test all(g -> length(components(g)) == 1,fixed)
    @test Set(vcat(fixed...)) == Set(sinks)
    @test sum(length,fixed) == length(sinks)
    @test C(4,10) in raw[1] && C(4,10) in fixed[2]
end
