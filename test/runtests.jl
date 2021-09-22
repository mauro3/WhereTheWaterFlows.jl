using WhereTheWaterFlows
const WWF = WhereTheWaterFlows
using Test

# test examples
module Test_Examples # use a module to avoid name-space pollution
plotyes = false
for fl in readdir("../examples/")
    include(joinpath("../examples/", fl))
end
end

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

function dem_one_point()
    dem = zeros(4,5) .* NaN
    dem[2,2] = 1.0
    coords = range(0, 1.0, length=size(dem,1)), range(0, 1.0, length=size(dem,2))
    return coords, dem
end

function dem_two_points()
    coords, dem = dem_one_point()
    dem[3,4] = 1.0
    return coords, dem
end

function dem_patho1()
    dem = [NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN 224973.93798828125 NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN 1.5264895434570312e6 780071.6101074219 611466.6186523438 68967.85522460938 NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN 340947.7197265625 1.7224892822265625e6 1.503103447265625e6 646871.5270996094 NaN NaN NaN NaN; NaN NaN NaN NaN NaN 493601.3024902344 NaN 953505.1794433594 1.33096515625e6 251280.46264648438 689901.5191650391 NaN NaN NaN NaN; NaN NaN NaN 1.210690463256836e6 2.2551471240234375e6 1.620159921875e6 1.9493891125488281e6 2.07443427734375e6 2.1659396948242188e6 10309.342041015625 NaN NaN NaN NaN NaN; NaN NaN NaN 1.4056021618652344e6 2.9576955078125e6 2.6524201977539062e6 2.3390813134765625e6 2.8080285302734375e6 2.6852572412109375e6 1.7304499780273438e6 NaN NaN NaN NaN NaN; NaN NaN NaN NaN 3.1962603833007812e6 3.308844228515625e6 3.3650431396484375e6 3.518876923828125e6 2.88937056640625e6 2.29438943359375e6 2.190683447265625e6 1.920630322265625e6 NaN NaN NaN; NaN NaN NaN NaN 2.521603028564453e6 3.1732824658203125e6 2.8601568969726562e6 3.4935072900390625e6 3.284232265625e6 2.7940991479492188e6 2.621620751953125e6 1.9656544458007812e6 NaN NaN NaN; NaN NaN NaN NaN 2.1466443212890625e6 2.5994008081054688e6 997498.37890625 3.0142672705078125e6 3.1849692797851562e6 2.6993770849609375e6 2.7337567944335938e6 2.1765316455078125e6 NaN NaN NaN; NaN NaN NaN NaN 1.4447062829589844e6 1.6352845727539062e6 297836.35498046875 2.5211809497070312e6 2.6398626635742188e6 2.2131605932617188e6 1.3856159045410156e6 NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN 932178.7219238281 NaN NaN NaN NaN NaN NaN; NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]
    coords = range(0, 1.0, length=size(dem,1)), range(0, 1.0, length=size(dem,2))
    return coords, dem
end

function dem_patho2()
    dem = [2.4669351831054688e6 2.4912734814453125e6 2.4511637768554688e6 2.297378525390625e6 2.1467358154296875e6 2.0279176049804688e6 1.8322703247070312e6 1.6419496728515625e6 1.3741556689453125e6 985379.6069335938 170602.13500976562; 2.6378740991210938e6 2.621620751953125e6 2.4854557543945312e6 2.4079890673828125e6 2.2196343212890625e6 1.9963872729492188e6 1.8812639379882812e6 1.6405933081054688e6 1.3537354772949219e6 937561.1877441406 NaN; 2.7794645190429688e6 2.6392191967773438e6 2.6666443579101562e6 2.4175318090820312e6 2.2348992724609375e6 2.0938221411132812e6 1.7992342114257812e6 1.6129501953125e6 1.1581694873046875e6 650434.7094726562 NaN; 2.7840773608398438e6 2.764429228515625e6 2.6759121826171875e6 2.5908404150390625e6 2.393361572265625e6 2.2053372387695312e6 2.128048779296875e6 1.8252124182128906e6 1.2706359448242188e6 820004.8461914062 NaN; 2.818317265625e6 2.876910830078125e6 2.7989625903320312e6 2.6884156298828125e6 2.3920689697265625e6 2.2358143334960938e6 2.211969111328125e6 1.9571274853515625e6 1.3319587353515625e6 864203.271484375 NaN; 2.7831515209960938e6 2.7405954125976562e6 2.728760107421875e6 2.5505696166992188e6 2.4580276049804688e6 2.0936141577148438e6 2.1094039184570312e6 1.9852256274414062e6 1.5260443835449219e6 868576.2139892578 NaN; 2.6441077563476562e6 2.666107099609375e6 2.6324500732421875e6 2.3455210571289062e6 2.2824994604492188e6 2.2010346118164062e6 2.149140135498047e6 1.762273427734375e6 1.3110776196289062e6 NaN NaN; 2.4332154907226562e6 2.3887434301757812e6 2.3717401684570312e6 2.1986926538085938e6 1.96104419921875e6 1.7805842163085938e6 1.8096322509765625e6 1.395735947265625e6 682117.8509521484 NaN NaN; 2.1785780810546875e6 2.1734258032226562e6 2.1239198754882812e6 1.974646943359375e6 1.6406056872558594e6 1.2984880700683594e6 1.258363583984375e6 726651.7492675781 NaN NaN NaN; 1.8946451831054688e6 1.81941337890625e6 1.7843751513671875e6 1.7247107556152344e6 1.2050103857421875e6 93141.61499023438 535295.9533691406 NaN NaN NaN NaN; 1.5536723583984375e6 1.3856159045410156e6 1.3403697265625e6 1.2850621618652344e6 793587.7758789062 NaN NaN NaN NaN NaN NaN]
    coords = range(0, 1.0, length=size(dem,1)), range(0, 1.0, length=size(dem,2))
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

@testset "diagonal gradient adjustment" begin
    dem = ones(4,4)
    dem[2,2] = 0.5
    dem[2,3] = 0.45
    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem);
    @test dir==[5  5  5  5
                5  8  1  5
                5  4  4  5 # without adjustment this row would be '5  7  4  5'
                5  5  5  5]
end


# Test high-level functions
@testset "DEM: dem1" begin
    dx = 0.9
    xs = -1.5:dx:1
    ys = -0.5:dx:3.0
    dem = dem1.(xs, ys')
    @test size(dem)==(length(xs), length(ys))

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, bnd_as_pits=false)
    @test area == [1.0 1.0 4.0 2.0; 1.0 1.0 1.0 1.0; 5.0 3.0 2.0 1.0]
    @test slen == [1 1 2 2; 1 1 1 1; 4 3 2 1]
    @test dir == Int8[5 8 5 5; 6 7 4 4; 5 2 2 2]
    @test nout == Bool[false true false false; true true true true; false true true true]
    @test nin == Int8[0 0 3 1; 0 0 0 0; 2 1 1 0]
    @test pits ==CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(3, 1), CartesianIndex(1, 3), CartesianIndex(1, 4)]

    # non default cellarea
    area, _ = WWF.waterflows(dem, fill(10.0, size(dem)), bnd_as_pits=false)
    @test area == [10.0 10.0 40.0 20.0; 10.0 10.0 10.0 10.0; 50.0 30.0 20.0 10.0]
end

@testset "DEM: peaks2" begin
    xs, dem = peaks2()
    ys = xs
    @test size(dem)==(length(xs), length(ys))
    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=false, bnd_as_pits=false);
    #plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 6
    @test maximum(slen)==86
    @test maximum(area)==4759
    @test length(unique(c))==6
    @test sum(c) == 33862
    @test sum(diff(c[:])) == 5
    @test sort(unique(c))[[1,end]] ==[1,6]
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=false)
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 4
    @test maximum(slen)==129
    @test maximum(area)==7702
    @test length(unique(c))==4
    @test sum(c) == 23855
    @test sum(diff(c[:])) == 3
    @test sort(unique(c))[[1,end]] ==[1,4]
    @test bnds isa Array{Array{CartesianIndex{2},1},1}
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color
end

@testset "DEM: peaks2_nan" begin
    xs, dem = peaks2_nan()
    ys = xs
    nanlocs = findall(isnan.(dem))
    @test size(dem)==(length(xs), length(ys))

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=false, bnd_as_pits=false);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 7
    @test maximum(slen)==87
    @test maximum(area[.!isnan.(area)])==4725
    @test length(unique(c))==8
    @test sum(c) == 36447
    @test sum(diff(c[:])) == 6
    @test sort(unique(c))[[1,end]] ==[0,7]
    @test all(c[nanlocs].==0)
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=false, bnd_as_pits=true);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 408
    @test maximum(slen)==66
    @test maximum(area[.!isnan.(area)])==4457
    @test length(unique(c))==409
    @test sum(c) == 2079380
    @test sum(diff(c[:])) == 407
    @test sort(unique(c))[[1,end]] ==[0,408]
    @test all(c[ [n for n in nanlocs if !(n in pits)]].==0)
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color

    # fill pits
    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=false);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 4
    @test maximum(slen)==142
    @test maximum(area[.!isnan.(area)])==7706
    @test length(unique(c))==5
    @test sum(c) == 23820
    @test sum(diff(c[:])) == 3
    @test sort(unique(c))[[1,end]] ==[0,4]
    @test bnds isa Array{Array{CartesianIndex{2},1},1}
    @test all(c[nanlocs].==0)
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=true);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 406
    @test maximum(slen)==75
    @test maximum(area[.!isnan.(area)])==4458
    @test length(unique(c))==407
    @test sum(c) == 2067155
    @test sum(diff(c[:])) == 405
    @test sort(unique(c))[[1,end]] ==[0,406]
    @test bnds isa Array{Array{CartesianIndex{2},1},1}
    @test all(c[ [n for n in nanlocs if !(n in pits)]].==0)
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color
end

@testset "DEM: peaks2_nan_edge" begin
    xs, dem = peaks2_nan_edge()
    ys = xs
    @test size(dem)==(length(xs), length(ys))

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=false, bnd_as_pits=false);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 6
    @test maximum(slen)==84
    @test maximum(area[.!isnan.(area)])==4658
    @test length(unique(c))==7
    @test sum(c) == 32486
    @test sum(diff(c[:])) == 0
    @test sort(unique(c))[[1,end]] ==[0,6]
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=false, bnd_as_pits=true);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 390
    @test maximum(slen)==64
    @test maximum(area[.!isnan.(area)])==4434
    @test length(unique(c))==391
    @test sum(c) == 2000894
    @test sum(diff(c[:])) == -1
    @test sort(unique(c))[[1,end]] ==[0,390]
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=false);
    #plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 1
    @test maximum(slen)==195
    @test maximum(area[.!isnan.(area)])==9605
    @test length(unique(c))==2
    @test sum(c) == 9605
    @test sum(diff(c[:])) == 0
    @test sort(unique(c))[[1,end]] ==[0,1]
    @test bnds isa Array{Array{CartesianIndex{2},1},1}
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=true);
    # plotarea_dem(xs, ys, dem, area, pits)
    @test length(pits) == 388
    @test maximum(slen)==81
    @test maximum(area[.!isnan.(area)])==4978
    @test length(unique(c))==389
    @test sum(c) == 1980928
    @test sum(diff(c[:])) == -1
    @test sort(unique(c))[[1,end]] ==[0,388]
    @test bnds isa Array{Array{CartesianIndex{2},1},1}
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color
end


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

# bug with non-empty set
@testset "non_empty" begin
    (xs, ys), dem = dem_one_point()
    mask = .!isnan.(dem)
    #should not error
    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=false)
    @test mask[pits[1]]
    @test length(pits)==1
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color
    #@test all(getindex.(Ref(mask), pits).==0) # tests that there are no interior pits left
    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=true)
    @test all(getindex.(Ref(mask), pits).==0)
    #area[isnan.(dem)] .= NaN; WWF.plotarea(xs, ys, area, pits)
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color

    (xs, ys), dem = dem_two_points()
    mask = .!isnan.(dem)
    #should not error
    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=false)
    @test mask[pits[1]]
    @test mask[pits[2]]
    @test length(pits)==2
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color
    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=true)
    @test all(getindex.(Ref(mask), pits).==0)
    #area[isnan.(dem)] .= NaN; WWF.plotarea(xs, ys, area, pits)
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color

    (xs, ys), dem = dem_patho1()
    mask = .!isnan.(dem)
    #should not error
    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=false)
    @test length(pits)==2
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color
    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=true)
    @test all(getindex.(Ref(mask), pits).==0)
    #area[isnan.(dem)] .= NaN; WWF.plotarea(xs, ys, area, pits)
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color

    # should not have pit in the interior
    (xs, ys), dem = dem_patho2()
    mask = .!isnan.(dem)
    # also mask border points
    mask[1,:] .= false; mask[end,:] .= false; mask[:,1] .= false; mask[:,end] .= false
    # should not have a pit in the interior
    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=false)
    @test all(getindex.(Ref(mask), pits).==0)
    @test all([c[pits[cc]]==cc  for cc=1:length(pits)]) # pit in catchment of same color
    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=true)
    @test all(getindex.(Ref(mask), pits).==0)
    #area[isnan.(dem)] .= NaN; WWF.plotarea(xs, ys, area, pits)
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
            @test catchment(dir, ij)==(c.==cc)
        end
    end
end

@testset "NOFLOWer" begin
    xs, dem = peaks2_nan_edge()
    ys = xs

    area, slen, dir, nout, nin, pits, c, bnds = WWF.waterflows(dem, drain_pits=false, bnd_as_pits=true);
    @test sum(c.==0) == 8
    # add a NOFLOWer cell
    # (Note: this will not be consistent with the routing
    # (i.e. flow is still into NOFLOWer, i.e. NOFLOWer should be a pit) and thus
    # leads to bongous results.  Well, let's just test the bongous results anyway,
    # to make sure they stay the same bongous.)
    ind = CartesianIndex(43,30)
    dir[ind] = WWF.NOFLOWer
    # and drain the pits
    WWF.drainpits!(dir, nin, nout, pits, c, bnds, dem)
    area, slen, c = WWF.flowrouting_catchments(dir, pits, ones(size(dem)))
    @test sum(c.==0) == 42 # these are bongous...

    area_flow, dir_flow = WWF.waterflows(dem, drain_pits=true, bnd_as_pits=true)[[1,3]];
    @test sum(area.!=area_flow) == 79
    @test sum(dir.!=dir_flow) == 1
end
