# AI disclosure: developed with an AI coding assistant; requires maintainer review.
module BlockSparseTests
using Test, Random, Serialization, WhereTheWaterFlows
const WWF = WhereTheWaterFlows
const WWFS = WWF.Subglacially
const WWFR = WWF.Randomly

occupied(A) = count(B -> B isa Array, WWFR.blocks(A))

@testset "Block-sparse catchment storage" begin
    sz = (123,117)
    counts = WWFR._init_catchments(sz,4,:blocksparse)
    dense = WWFR._init_catchments(sz,4,:dense)
    @test counts isa AbstractArray{Float16,3}
    @test size(counts)==(sz...,4)
    @test axes(counts)==axes(dense)
    @test size(WWFR.blocks(counts))==(3,3,4)
    @test occupied(counts)==0
    @test WWFR.blocks(counts)[1,1,1] === WWFR.blocks(counts)[2,2,4]
    @test counts[123,117,4]==0
    @test counts[CartesianIndex(51,51,2)]==0
    @test Array(counts)==dense
    @test_throws ArgumentError WWFR._init_catchments(sz,1,:unknown)
    @test_throws DimensionMismatch WWFR._accumulate_catchment!(counts,falses(2,3),1)
    @test_throws BoundsError WWFR._accumulate_catchment!(counts,falses(sz),5)
    @test WWFR._accumulate_catchment!(counts,falses(sz),1) === counts
    @test occupied(counts)==0

    mask = falses(sz)
    for I in ((1,1),(50,50),(51,51),(100,100),(101,101),(123,117))
        mask[I...] = true
    end
    for draw in 1:3
        if draw==2
            # New support in a formerly empty tile; earlier memberships must persist.
            mask .= false
            mask[51,1] = true
            mask[123,117] = true
        end
        # Duplicate groups overlap but must own independent dense blocks.
        for group in (1,2,3)
            selection = group==2 ? Matrix(mask) : mask
            WWFR._accumulate_catchment!(counts,selection,group)
            WWFR._accumulate_catchment!(dense,selection,group)
        end
        @test counts==dense
        @test Array(counts[:,:,4])==zeros(Float16,sz)
        @test all(B -> B isa WWFR.Zeros,view(WWFR.blocks(counts),:,:,4))
    end
    @test counts[1,1,1]==1
    @test counts[51,1,1]==2
    @test counts[123,117,1]==3
    @test size(WWFR.blocks(counts)[3,3,1])==(23,17,1)
    @test WWFR.blocks(counts)[1,1,1] !== WWFR.blocks(counts)[1,1,3]
    first_tile = WWFR.blocks(counts)[2,1,1]
    # Warmed updates of existing support must not allocate new tile storage.
    WWFR._accumulate_catchment!(counts,mask,1)
    WWFR._accumulate_catchment!(dense,mask,1)
    bytes = @allocated WWFR._accumulate_catchment!(counts,mask,1)
    WWFR._accumulate_catchment!(dense,mask,1)
    @test bytes<1024
    @test WWFR.blocks(counts)[2,1,1] === first_tile
    before = occupied(counts)
    @test WWFR._normalize_catchments!(counts,5) === counts
    WWFR._normalize_catchments!(dense,5)
    @test counts==dense
    @test occupied(counts)==before
    @test (@allocated WWFR._normalize_catchments!(counts,1))<1024
    @test WWFR.blocks(counts)[2,1,1] === first_tile
    slice = view(counts,20:3:110,10:4:110,1)
    @test Array(slice)==dense[20:3:110,10:4:110,1]
    @test Array(copy(counts))==dense
    @test Array(counts .* 2)==dense .* 2
    io=IOBuffer()
    serialize(io,counts)
    seekstart(io)
    loaded=deserialize(io)
    @test loaded==counts
    @test occupied(loaded)==before
    @test all(B -> B isa WWFR.Zeros,view(WWFR.blocks(loaded),:,:,4))

    empty = WWFR._init_catchments(sz,0,:blocksparse)
    @test size(empty)==(sz...,0)
    @test isempty(empty)
    @test size(Array(empty))==(sz...,0)
    @test WWFR._normalize_catchments!(empty,3) === empty
    exact = WWFR._init_catchments((2,3),1,:blocksparse)
    for _ in 1:2048
        WWFR._accumulate_catchment!(exact,trues(2,3),1)
    end
    @test all(==(Float16(2048)),exact)
    WWFR._normalize_catchments!(exact,2048)
    @test all(isone,exact)
    @test_throws ErrorException WWFR.map_mc(nothing,nothing,nothing,2049;progressmeter=false)
end

@testset "Native-shaped catchments without dense values" begin
    counts=WWFR._init_catchments((11061,9001),198,:blocksparse)
    empty_bytes=Base.summarysize(counts)
    @test size(counts)==(11061,9001,198)
    @test size(WWFR.blocks(counts))==(222,181,198)
    @test empty_bytes<70*1024^2
    mask=falses(11061,9001)
    mask[51,101]=true
    WWFR._accumulate_catchment!(counts,mask,1)
    mask[51,101]=false
    mask[11061,9001]=true
    WWFR._accumulate_catchment!(counts,mask,198)
    @test counts[51,101,1]==1
    @test counts[11061,9001,198]==1
    @test counts[51,101,198]==0
    @test occupied(counts)==2
    occupied_bytes=Base.summarysize(counts)
    @test occupied_bytes-empty_bytes<20_000
    @test (@allocated WWFR._normalize_catchments!(counts,2))<1024
    @test counts[51,101,1]==0.5
    @test counts[11061,9001,198]==0.5
    @test occupied(counts)==2
    @info "Native-shaped catchment storage (object bytes, not RSS)" empty_bytes occupied_bytes
end

@testset "Dense/block-sparse reducer equivalence" begin
    nx,ny=63,57
    bed=[1000.0-10i+0.5j for i in 1:nx,j in 1:ny]
    surface=bed .+ 100
    source=fill(1e-8,nx,ny)
    ff=fill(0.97,nx,ny)
    fixed=WWFR.Uncertainty()
    uncertain=WWFR.Uncertainty(absuc=1.0,correlation_length=400.0)
    baseline=WWF.waterflows(bed)
    groups=[baseline.sinks,[CartesianIndex(51,51)],baseline.sinks,CartesianIndex{2}[]]
    for kind in (:subaerial,:subglacial), gamma in (0.0,-0.31)
        build(storage) = kind==:subaerial ?
            WWFR.make_fns_subaerial(200.0,bed,uncertain,source,fixed,groups;
                catchment_storage=storage) :
            WWFR.make_fns_subglacial(200.0,surface,fixed,bed,uncertain,ff,fixed,
                source,fixed,groups;gamma,catchment_storage=storage)
        model,sample,reduce_sparse! = build(:blocksparse)
        _,_,reduce_dense! = build(:dense)
        sparse,dense=reduce_sparse!(),reduce_dense!()
        Random.seed!(731)
        for _ in 1:3
            res=model(sample()...)
            if kind==:subglacial
                @test all(M -> M isa BitMatrix,res[2].sink_catchments.masks)
            end
            reduce_sparse!(sparse,res)
            reduce_dense!(dense,res)
            @test Array(sparse.catchments)==dense.catchments
        end
        parent_before=WWFR.blocks(sparse.catchments)[1,1,1]
        reduce_sparse!(sparse)
        reduce_dense!(dense)
        @test WWFR.blocks(sparse.catchments)[1,1,1] === parent_before
        for field in keys(dense)
            actual,expected=getproperty(sparse,field),getproperty(dense,field)
            if field==:n_samples
                @test actual[]==expected[]==3
            elseif field==:catchments
                @test Array(actual)==expected
            else
                @test isequal(actual,expected)
            end
        end
        @test all(B -> B isa WWFR.Zeros,view(WWFR.blocks(sparse.catchments),:,:,4))
    end
end
end
