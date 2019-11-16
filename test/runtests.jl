using WhereTheWaterFlows
const WWF = WhereTheWaterFlows
using Test


"An artificial DEM"
function dem1(x, y; withpit=false, randfac=0.0)
    out = - (x^2 - 1)^2 - (x^2*y - x - 1)^2 + 6 + 0.1*x + 3*y
    if withpit
        out -= 2*exp(-(x^2 + y^2)*50)
    end
    out += randfac*randn()
    return out<0 ? 0.0 : out
end

const _offset = [0.2, 0.5]

"DEM with a few maxs and mins"
function peaks(n=100)
    coords = range(-pi+_offset[1], pi-_offset[2], length=n)
    return coords, sin.(coords) .* cos.(coords')
end

"DEM with a few more maxs and mins"
function peaks2(n=100, randfac=0.0)
    coords = range(-pi+_offset[1], pi-_offset[2], length=n)
    return coords, sin.(coords) .* cos.(coords') .-
        0.7*(sin.(coords.+1) .* cos.(coords')).^8 .+
        0.01 * (coords.+coords') .+
        randfac*randn(n,n) # 0.02
end

"DEM with a few more maxs and mins.  NaN masked middle."
function peaks2_nan(n=100, randfac=0.0)
    coords = range(-pi+_offset[1], pi-_offset[2], length=n)
    dem = sin.(coords) .* cos.(coords') .-
        0.7*(sin.(coords.+1) .* cos.(coords')).^8 .+
        0.01 * coords.*coords' .+
        randfac*randn(n,n) # 0.02
    dem[(n÷2:n÷2+2).-8, n÷2+1:n÷2+4] .= NaN
    return coords, dem
end

"DEM with a few more maxs and mins.  NaN masked edge."
function peaks2_nan_edge(n=100, randfac=0.0)
    coords = range(-pi+_offset[1], pi-_offset[2], length=n)
    dem = sin.(coords) .* cos.(coords') .-
        0.7*(sin.(coords.+1) .* cos.(coords')).^8 .+
        0.01 * (coords.+coords') .+
        randfac*randn(n,n) # 0.02
    outflow = dem[n÷2, 1]
    dem[1, :] .= NaN
    dem[end, :] .= NaN
    dem[:, 1] .= NaN
    dem[:, end] .= NaN
    dem[n÷2, 1] = outflow
    return coords, dem
end

# Test low-level fns

# flowsinto
# """
# Tests whether a cell `J` with flowdir `dirJ` flows into cell `I`.
# """
# flowsinto(J::CartesianIndex, dirJ, I::CartesianIndex) = ind2dir(I-J) == dirJ
@testset "flowsinto" begin
    J = CartesianIndex(5,6)
    dirs = 0:9
    for i=[3,7], j=4:8, dir in dirs
        @test_throws BoundsError WWF.flowsinto(J, dir, CartesianIndex(i,j))
    end
    for i=3:7, j=[4,8], dir in dirs
        @test_throws BoundsError WWF.flowsinto(J, dir, CartesianIndex(i,j))
    end
    trueflow = Dict(1 => J+CartesianIndex(-1,-1),
                    2 => J+CartesianIndex(0,-1),
                    3 => J+CartesianIndex(1,-1),
                    4 => J+CartesianIndex(-1,0),
                    5 => J+CartesianIndex(0,0),
                    6 => J+CartesianIndex(1,0),
                    7 => J+CartesianIndex(-1,1),
                    8 => J+CartesianIndex(0,1),
                    9 => J+CartesianIndex(1,1)
                    )
    for (td, I) in trueflow
        for d in dirs
            if d==td
                @test WWF.flowsinto(J, d, I)
            else
                @test !WWF.flowsinto(J, d, I)
            end
        end
    end
end

@test WWF.dirnums==reverse([ 7 8 9
                             4 5 6
                             1 2 3]', dims=2)
tmp = convert(Matrix, WWF.dirnums)

for i=1:9
    dem = ones(3,3)
    dem[i] = 0
    dir, nout, nin, pits = WWF.d8dir_feature(dem, false)

    @test dir[2,2] == WWF.dirnums[i]
    @test dem[WWF.dir2ind(dir[2,2])+CartesianIndex(2,2)]==0

    @test WWF.d8dir_feature(dem, false)[1][2,2] == i
    @test maximum(nout)<=1
    @test minimum(nout)>=0
    @test maximum(nin)<=8
    @test minimum(nin)>=0
    @test sum(nout)==sum(nin)

    @test sum(dem[WWF.iterate_D9(CartesianIndex(2,2), dem)]) == 8

end

@testset "dir2ind" begin
    @test WWF.dir2ind(1) == CartesianIndex(-1,-1)
    @test WWF.dir2ind(2) == CartesianIndex(0,-1)
    @test WWF.dir2ind(3) == CartesianIndex(1,-1)
    @test WWF.dir2ind(4) == CartesianIndex(-1,0)
    @test WWF.dir2ind(5) == CartesianIndex(0,0)
    @test WWF.dir2ind(6) == CartesianIndex(1,0)
    @test WWF.dir2ind(7) == CartesianIndex(-1,1)
    @test WWF.dir2ind(8) == CartesianIndex(0,1)
    @test WWF.dir2ind(9) == CartesianIndex(1,1)
    for i=1:9
        @test i==WWF.ind2dir(WWF.dir2ind(i))
    end
end


# Test high-level functions
@testset "DEM: dem1" begin
    dx = 0.9
    xs = -1.5:dx:1
    ys = -0.5:dx:3.0
    dem = dem1.(xs, ys')
    @test size(dem)==(length(xs), length(ys))

    area, slen, dir, nout, nin, pits = WWF.waterflows(dem, bnd_as_pits=false)
    @test area == [1.0 1.0 5.0 1.0; 1.0 1.0 1.0 1.0; 5.0 3.0 2.0 1.0]
    @test slen == [1 1 2 1; 1 1 1 1; 4 3 2 1]
    @test dir == Int8[5 8 5 5; 6 7 4 1; 5 2 2 2]
    @test nout == Bool[false true false false; true true true true; false true true true]
    @test nin == Int8[0 0 4 0; 0 0 0 0; 2 1 1 0]
    @test pits ==CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(3, 1), CartesianIndex(1, 3), CartesianIndex(1, 4)]
end

@testset "DEM: peaks2" begin
    xs, dem = peaks2()
    ys = xs
    @test size(dem)==(length(xs), length(ys))
    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=false);
    #plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 6
    @test maximum(slen)==67
    @test maximum(area)==4831
    @test length(unique(c))==6
    @test sum(c) == 33876
    @test sum(diff(c[:])) == 5
    @test sort(unique(c))[[1,end]] ==[1,6]

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 4
    @test maximum(slen)==118
    @test maximum(area)==7714
    @test length(unique(c))==4
    @test sum(c) == 21990
    @test sum(diff(c[:])) == 2
    @test sort(unique(c))[[1,end]] ==[1,4]
    @test bnds isa Array{Array{CartesianIndex,1},1}
end

@testset "DEM: peaks2_nan" begin
    xs, dem = peaks2_nan()
    ys = xs
    @test size(dem)==(length(xs), length(ys))

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=false, bnd_as_pits=false);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 7
    @test maximum(slen)==67
    @test maximum(area[.!isnan.(area)])==4625
    @test length(unique(c))==8
    @test sum(c) == 36607
    @test sum(diff(c[:])) == 6
    @test sort(unique(c))[[1,end]] ==[0,7]

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=false, bnd_as_pits=true);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 408
    @test maximum(slen)==66
    @test maximum(area[.!isnan.(area)])==4114
    @test length(unique(c))==409
    @test sum(c) == 2083959
    @test sum(diff(c[:])) == 407
    @test sort(unique(c))[[1,end]] ==[0,408]

    # fill pits
    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=false);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 4
    @test maximum(slen)==116
    @test maximum(area[.!isnan.(area)])==7707
    @test length(unique(c))==5
    @test sum(c) == 24117
    @test sum(diff(c[:])) == 3
    @test sort(unique(c))[[1,end]] ==[1,5]
    @test bnds isa Array{Array{CartesianIndex,1},1}

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=true);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 406
    @test maximum(slen)==75
    @test maximum(area[.!isnan.(area)])==4115
    @test length(unique(c))==407
    @test sum(c) == 1444498
    @test sum(diff(c[:])) == 406
    @test sort(unique(c))[[1,end]] ==[1,407]
    @test bnds isa Array{Array{CartesianIndex,1},1}
end

@testset "DEM: peaks2_nan_edge" begin
    xs, dem = peaks2_nan_edge()
    ys = xs
    @test size(dem)==(length(xs), length(ys))

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=false, bnd_as_pits=false);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 6
    @test maximum(slen)==66
    @test maximum(area[.!isnan.(area)])==4727
    @test length(unique(c))==7
    @test sum(c) == 32498
    @test sum(diff(c[:])) == 0
    @test sort(unique(c))[[1,end]] ==[0,6]

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=false, bnd_as_pits=true);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 390
    @test maximum(slen)==65
    @test maximum(area[.!isnan.(area)])==4372
    @test length(unique(c))==391
    @test sum(c) == 2003347
    @test sum(diff(c[:])) == -1
    @test sort(unique(c))[[1,end]] ==[0,390]


    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=false);
    #plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 1
    @test maximum(slen)==165
    @test maximum(area[.!isnan.(area)])==9605
    @test length(unique(c))==2
    @test sum(c) == 19605
    @test sum(diff(c[:])) == 0
    @test sort(unique(c))[[1,end]] ==[1,2]
    @test bnds isa Array{Array{CartesianIndex,1},1}

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=true);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 388
    @test maximum(slen)==81
    @test maximum(area[.!isnan.(area)])==4831
    @test length(unique(c))==389
    @test sum(c) == 1372917
    @test sum(diff(c[:])) == 98
    @test sort(unique(c))[[1,end]] ==[1,389]
    @test bnds isa Array{Array{CartesianIndex,1},1}
end


@testset "fill_dem" begin
    # TODO: need more fill tests
    dx = 0.1
    xs = -1.5:dx:1
    ys = -0.5:dx:3.0
    dem = dem1.(xs, ys', withpit=true)
    area, slen, dir, nout, nin, pits = WWF.waterflows(dem)
    demf = WWF.fill_dem(dem, pits, dir)
    @test sum(demf.-dem) ≈ 2.1499674517313414
    @test sum(demf.-dem .> 0) == 5
end
