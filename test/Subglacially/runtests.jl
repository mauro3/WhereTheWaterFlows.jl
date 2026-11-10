module Subglacially
using WhereTheWaterFlows
const WWF = WhereTheWaterFlows
const WWFS = WWF.Subglacially
using Test

include("catchments.jl")

function test_against_reference(relpath::AbstractString, value)
    path = joinpath(@__DIR__, relpath)
    isfile(path) || error("Missing reference file: $relpath")
    @test strip(read(path, String)) == strip(string(value))
end

const _offset = [0.2, 0.5]
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

"Smooth ice-surface DEM"
function icesurf_nan_edge(n=100, randfac=0.0)
    coords = range(-pi+_offset[1], pi-_offset[2], length=n)
    f = 0.01
    #dem = -((f*coords/pi).^2 .+ ((f*coords/pi).^2)') .+ 2
    dem = 0.05*coords .+ (0.05*coords') .+ 4

    outflow = dem[n÷2, 1]
    dem[1, :] .= NaN
    dem[end, :] .= NaN
    dem[:, 1] .= NaN
    dem[:, end] .= NaN
    dem[n÷2, 1] = outflow
    return coords, dem
end


# Just some random indices for reference tests
inds =[644, 961, 1281, 1295, 1486, 1669, 1881, 1884, 2250, 2466, 2832, 2881, 3129, 3488, 3531,
       4103, 4171, 4246, 4407, 4763, 4839, 4907, 5568, 5945, 6318, 6691, 6945, 7278, 7602, 7673,
       7767, 8822, 8864, 8951, 9543, 9709, 9716, 9850, 9898, 9992]

function next_to_outer_boundary(ar, I::CartesianIndex)
    i,j = I.I
    iend, jend = size(ar)
    if i==2 || i==iend-1 || j==2 || j==jend-1
        return true
    else
        return false
    end
end


@testset "Subglacial" begin
    n = 100
    xs, beddem = peaks2_nan_edge(n)
    xs, surfdem = icesurf_nan_edge(n)
    xs = xs.*10000 # to make a decent size ice sheet
    dx = step(xs)
    beddem *= 20
    surfdem *= 100
    ys = xs
    thick = (surfdem .- beddem)
    phi = thick * 910 * 9.8 + beddem * 1000 * 9.8
    phim =  beddem * 1000 * 9.8
    floatfrac = 1
    source = fill!(similar(surfdem), 1)/dx^2

    for dp=[true,false], bap=[true,false], gamma=[0, WWFS.GAMMA], avoid_sc=[true,false]
        # dp=[true,false][1]; bap=[true,false][1]; gamma=[0, WWFS.GAMMA][1]; avoid_sc=[true,false][2]
        # @show dp, bap, gamma, avoid_sc
        out = WWFS.waterflows_subglacial(surfdem, beddem, dx, floatfrac, source; drain_pits=dp, bnd_as_sink=true, nan_as_sink=bap, gamma=gamma, avoid_sc=avoid_sc, rhoi=910)
        (;slen, dir, nout, nin, sinks, pits, c, bnds, area, phi) = out.routing
        area = area.total
        (;sc_locs, kappas, dir_og) = out.pressmelt

        for (name, var) in zip((:area, :slen, :dir, :nout, :nin, :c, :sc_locs, :kappas, :dir_og), (area, slen, dir, nout, nin, c, sc_locs, kappas, dir_og))
            # @show name
            test_against_reference("reftest-files/wwfs-$gamma-$(avoid_sc)-$dp-$bap-$name.test", var[inds])
        end
        @test length(bnds)==length(pits)
        inds_ = round.(Int, range(1,length(sinks),length=20))
        test_against_reference("reftest-files/wwfs-$gamma-$(avoid_sc)-$dp-$bap-sinks.test", sinks[inds_])

        if gamma==0
            area_, slen_, dir_, nout_, nin_, sinks_, pits_, c_, bnds_ = WWF.waterflows(phi, drain_pits=dp, bnd_as_sink=true, nan_as_sink=bap);
            # not true because of extra melt: @test area==area_
            @test slen==slen_
            @test dir==dir_
            @test nout==nout_
            @test nin==nin_
            @test sinks==sinks_
            @test pits==pits_
            @test c==c_
            @test bnds==bnds_
        end
    end
end

@testset "Automatic ice-free feedback gate retains prescribed source" begin
    WWFR = WWF.Randomly
    bed = [2000. - 1000i + 5j for i in 1:7,j in 1:9]
    surface = bed .+ 50.
    source = fill(1e-8,7,9)
    original_source = copy(source)
    plain = WWF.waterflows(bed,source .* 500^2)
    for gamma in (0.,-0.31)
        noice = WWFS.waterflows_subglacial(bed,bed,500.,1.,source;gamma)
        @test noice.routing.dir == plain.dir
        @test noice.routing.area.total == plain.area
        @test all(iszero,noice.routing.area.dissipation_melt_rate)
        @test all(iszero,noice.routing.area.pressure_melt_rate)
        @test sum(noice.routing.area.total[I] for I in noice.routing.sinks) ≈ sum(source)*500^2
        @test source == original_source
        ice = WWFS.waterflows_subglacial(surface,bed,500.,1.,source;gamma)
        @test any(!iszero,ice.routing.area.dissipation_melt_rate)
    end
    fixed = WWFR.Uncertainty()
    (;model,initialize,aggregate!,finalize!) = WWFR.make_fns_subglacial(500.,
        surface,fixed,bed,fixed,ones(7,9),fixed,source,fixed,[];gamma=-.31)
    # Factory reference has ice, but the current model call does not.
    fields,noice = model(bed,bed,ones(7,9),source)
    @test fields.source === source
    @test all(iszero,noice.routing.area.dissipation_melt_rate)
    current_sample() = (bed,bed,ones(7,9),source)
    aggr = WWFR.map_mc(model,current_sample,initialize,aggregate!,finalize!,2;progressmeter=false)
    @test all(iszero,aggr.melt_rate)
    @test aggr.areas_total ≈ plain.area rtol=1e-6
    @test aggr.n_samples[]==2
end

@testset "Surface smoothing half-window stays within thickness" begin
    grid = 0.:500.:3000.
    surface = [1000. + i^2 + j^2 for i in 1:7, j in 1:7]
    for H in (250.,300.,499.999,500.,750.,999.999,1000.)
        bed = surface .- H
        result = WWFS.smooth_surface(grid,grid,surface,bed,1.)
        radius = floor(Int,H/500.)
        @test result[4,4] ≈ sum(surface[4-radius:4+radius,4-radius:4+radius]) / (2radius+1)^2
        @test radius * 500 <= H
        H < 500 && @test result == surface
    end
    @test WWFS.smooth_surface(grid,grid,surface,surface .- 900.,0.5) == surface
end

@testset "One-based input grids" begin
    WWFS = WWF.Subglacially
    dem = [3.0 2 1; 4 5 2; 6 7 3]
    shifted = WWF.Origin(0, 0)(copy(dem))
    source = ones(3, 3)
    for args in ((shifted, dem, 1.0), (dem, shifted, 1.0),
                 (dem, dem, 1.0, shifted), (dem, dem, 1.0, 1, shifted),
                 (dem, dem, 1.0, 1, source, shifted .> 0))
        @test_throws ArgumentError WWFS.waterflows_subglacial(args...)
    end
    @test_throws ArgumentError WWFS.boxcar(shifted, ones(Int, 3, 3))
end

end
