using WhereTheWaterFlows
using Test
const WWF = WhereTheWaterFlows

"DEM with a few maxs and mins"
function peaks(n=100)
    coords = range(-pi, pi, length=n)
    return coords, sin.(coords) .* cos.(coords')
end

"DEM with a few more maxs and mins"
function peaks2(n=100, randfac=0.0)
    coords = range(-pi, pi, length=n)
    return coords, sin.(coords) .* cos.(coords') .-
        0.7*(sin.(coords.+1) .* cos.(coords')).^8 .+
        randfac*randn(n,n) # 0.02
end

"Another artificial DEM"
function ele(x, y; withpit=false, randfac=0.0)
    out = - (x^2 - 1)^2 - (x^2*y - x - 1)^2 + 6 + 0.1*x + 3*y
    if withpit
        out -= 2*exp(-(x^2 + y^2)*50)
    end
    out += randfac*randn(size(out))
    return out<0 ? 0.0 : out
end

# Test low-level fns

@test WWF.dirnums==reverse([ 7 8 9
                            4 5 6
                            1 2 3]', dims=2)
tmp = convert(Matrix, WWF.dirnums)
tmp[2,2] = 0
@test WWF.dirnums_0==tmp

for i=1:9
    dem = ones(3,3)
    dem[i] = 0
    dir, nout, nin, pits = WWF.d8dir_feature(dem)

    @test dir[2,2] == WWF.dirnums[i]
    @test dem[WWF.dir2ind(dir[2,2])+CartesianIndex(2,2)]==0

    @test WWF.d8dir_feature(dem)[1][2,2] == i
    @test maximum(nout)<=1
    @test minimum(nout)>=0
    @test maximum(nin)<=8
    @test minimum(nin)>=0
    @test sum(nout)==sum(nin)

    @test sum(dem[WWF.iterate_D9(CartesianIndex(2,2), dem)]) == 8

end

@test WWF.dir2ind(1) == CartesianIndex(-1,-1)
@test WWF.dir2ind(2) == CartesianIndex(0,-1)
@test WWF.dir2ind(3) == CartesianIndex(1,-1)
@test WWF.dir2ind(4) == CartesianIndex(-1,0)
@test WWF.dir2ind(5) == CartesianIndex(0,0)
@test WWF.dir2ind(6) == CartesianIndex(1,0)
@test WWF.dir2ind(7) == CartesianIndex(-1,1)
@test WWF.dir2ind(8) == CartesianIndex(0,1)
@test WWF.dir2ind(9) == CartesianIndex(1,1)


# Test high-level functions
dx = 0.9
xs = -1.5:dx:1
ys = -0.5:dx:3.0
dem = ele.(xs, ys')
@test size(dem)==(length(xs), length(ys))


area, slen, dir, nout, nin, pits = WWF.waterflows(dem)
@test area == [1.0 1.0 5.0 1.0; 1.0 1.0 1.0 1.0; 5.0 3.0 2.0 1.0]
@test slen == [1 1 2 1; 1 1 1 1; 4 3 2 1]
@test dir == Int8[5 8 5 5; 6 7 4 1; 5 2 2 2]
@test nout == Bool[false true false false; true true true true; false true true true]
@test nin == Int8[0 0 4 0; 0 0 0 0; 2 1 1 0]
@test pits ==CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(3, 1), CartesianIndex(1, 3), CartesianIndex(1, 4)]

# Tests with peaks2
xs, dem = peaks2()
ys = xs
@test size(dem)==(length(xs), length(ys))
area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, fillpits=false);
@test length(pits) == 12
@test maximum(slen)==61
@test maximum(area)==3707
@test length(unique(c))==12
@test sum(c) == 60221
@test sum(diff(c[:])) == 11
@test sort(unique(c))[[1,end]] ==[1,12]
area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, fillpits=true);
@test length(pits) == 8
@test maximum(slen)==113
@test maximum(area)==6492
@test length(unique(c))==8
@test sum(c) == 38590
@test sum(diff(c[:])) == 5
@test sort(unique(c))[[1,end]] ==[1,8]
@test bnds isa Array{Array{CartesianIndex,1},1}
#WWF.plotarea(xs,ys,area,pits)


# check DEM fill
dx = 0.1
xs = -1.5:dx:1
ys = -0.5:dx:3.0
dem = ele.(xs, ys', withpit=true)
area, slen, dir, nout, nin, pits = WWF.waterflows(dem)
demf = WWF.fill_dem(dem, pits, dir)
@test sum(demf.-dem) ≈ 2.1499674517313414
@test sum(demf.-dem .> 0) == 5


demf = WWF.fill_dem(dem, pits, dir)
WWF.heatmap(xs,ys,demf)
