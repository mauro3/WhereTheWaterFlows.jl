module PublicAPI

using Test
using WhereTheWaterFlows

const WWF = WhereTheWaterFlows
const WWFS = WWF.Subglacially
const WWFR = WWF.Randomly

@testset "Root public API" begin
    @test Set(filter(name -> Base.isexported(WWF, name), names(WWF))) == Set((
        :WhereTheWaterFlows, :waterflows, :fill_dem, :catchment, :catchment!,
        :catchments, :catchment_flux, :prune_catchments,
        :plt_dir, :plt_dir!, :plt_area, :plt_area!, :plt_sinks, :plt_sinks!,
        :plt_bnds, :plt_bnds!, :plt_catchments, :plt_catchments!, :plt_it))
    for name in (:Subglacially, :Randomly, :PIT, :SINK, :BARRIER,
                 :dir2ind, :dir2vec, :ind2dir, :showme)
        @test Base.ispublic(WWF, name)
        @test !Base.isexported(WWF, name)
    end
    for name in (:d8dir_feature, :flowrouting_catchments, :make_flowfeatures,
                 :make_boundaries, :drainpits!, :dirnums, :cartesian)
        @test !Base.ispublic(WWF, name)
    end
end

@testset "Submodule public API" begin
    @test Set(filter(name -> Base.isexported(WWFR, name), names(WWFR))) == Set((:Randomly, :Uncertainty, :map_mc,
        :make_fns_subaerial, :make_fns_subglacial, :gaussian_kernel, :exponential_kernel))
    @test Set(filter(name -> Base.isexported(WWFS, name), names(WWFS))) == Set((:Subglacially, :waterflows_subglacial,
        :phi_fn, :lake_depth, :smooth_surface, :catchment_sinks))

    for name in (:mask_contiguous, :get_boundary_cells, :catchment_sink_components,
                 :connect_catchment_sinks!, :walk_gl, :RHOI, :RHOW, :GRAV, :LFUSION, :GAMMA)
        @test Base.ispublic(WWFS, name)
        @test !Base.isexported(WWFS, name)
    end
    for name in (:make_sampler, :make_field_realization, :make_grf_sampler,
                 :find_next_2357, :pad, :depad, :distance_squared, :make_kernel)
        @test !Base.ispublic(WWFR, name)
    end
    for name in (:d8dir_pressmelt, :make_feedback_fn, :melt_rates, :boxcar,
                 :signed_distance, :check_catchment_sinks, :walk_gl_inside, :walk_gl_outside)
        @test !Base.ispublic(WWFS, name)
    end
end

end
