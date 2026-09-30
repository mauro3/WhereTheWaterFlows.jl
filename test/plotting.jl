using Test, WhereTheWaterFlows
using CairoMakie

@testset "Reactive catchment filtering" begin
    original = [1 1 2; 1 1 3; 4 4 3]
    labels = Observable(copy(original))
    minsize = Observable(2)
    fp = plt_catchments(1:3, 1:3, labels; minsize)
    hm = only(fp.plot.plots)
    @test hm[3][] == [1 1 0; 1 1 2; 3 3 2]
    @test Tuple(hm.colorrange[]) == (1, 3)
    @test labels[] == original

    minsize[] = 3
    @test hm[3][] == [1 1 0; 1 1 0; 0 0 0]
    @test Tuple(hm.colorrange[]) == (1, 2)
    replacement = [0 5 5; 8 8 8; 8 0 5]
    labels[] = replacement
    @test hm[3][] == [0 1 1; 2 2 2; 2 0 1]
    @test labels[] == replacement

    minsize[] = 0
    @test hm[3][] == replacement
    @test Tuple(hm.colorrange[]) == (1, 8)
    minsize[] = 10
    @test hm[3][] == zeros(Int, 3, 3)
    @test Tuple(hm.colorrange[]) == (1, 2)
end

@testset "Boundary and sink point recipes" begin
    for (standalone, overlay) in ((plt_bnds, plt_bnds!), (plt_sinks, plt_sinks!))
        x = Observable([10.0, 20.0, 30.0])
        y = [40.0, 50.0, 60.0]
        indices = Observable([CartesianIndex(1, 2)])
        fp = standalone(x, y, indices)
        child = only(fp.plot.plots)
        @test child isa Makie.Scatter
        @test child[1][] == [Point2f(10, 50)]

        fig = Figure()
        ax = Axis(fig[1, 1])
        other = Figure()
        otherax = Axis(other[1, 1])
        p = overlay(ax, x, y, indices)
        @test p in ax.scene.plots
        @test isempty(otherax.scene.plots)
        overlay_child = only(p.plots)
        @test overlay_child[1][] == child[1][]

        # Change both coordinates and point count, including an empty list.
        indices[] = [CartesianIndex(2, 1), CartesianIndex(3, 3)]
        @test child[1][] == overlay_child[1][] == [Point2f(20, 40), Point2f(30, 60)]
        x[] = [11.0, 21.0, 31.0]
        @test child[1][] == overlay_child[1][] == [Point2f(21, 40), Point2f(31, 60)]
        indices[] = CartesianIndex{2}[]
        @test isempty(child[1][])
        @test isempty(overlay_child[1][])
    end
end

@testset "Area colorbar placement" begin
    glb = Makie.GridLayoutBase
    x = y = 1:3
    area = reshape(collect(1.0:9.0), 3, 3)
    fig = Figure()
    ax1 = Axis(fig[1, 1])
    ax2 = Axis(fig[1, 2])
    other = Figure()
    otherax = Axis(other[1, 1])
    p1 = plt_area!(ax1, x, y, area; colorbar_label="Area")
    p2 = plt_area!(ax2, x, y, area; colorbar_kwargs=(label="Custom",))
    bars = filter(b -> b isa Colorbar, fig.content)
    @test length(bars) == 2
    @test isempty(filter(b -> b isa Colorbar, other.content))
    @test Set(b.label[] for b in bars) == Set(["Area", "Custom"])
    @test p1 in ax1.scene.plots
    @test p2 in ax2.scene.plots
    for (ax, col) in ((ax1, 1), (ax2, 2))
        panel = glb.gridcontent(ax).parent
        @test glb.gridcontent(panel).parent === fig.layout
        @test glb.gridcontent(panel).span.cols == col:col
        @test count(b -> glb.gridcontent(b).parent === panel, bars) == 1
    end
    # Opting out preserves the original layout, and overlays still target the axis.
    plt_area!(otherax, x, y, area; colorbar=false)
    @test glb.gridcontent(otherax).parent === other.layout
    @test isempty(filter(b -> b isa Colorbar, other.content))
    overlay = scatter!(ax1, [1, 2], [2, 3])
    @test overlay in ax1.scene.plots

    standalone = plt_area(x, y, area)
    @test count(b -> b isa Colorbar, standalone.figure.content) == 1
    gridfig = Figure()
    ap = plt_area(gridfig[2, 1], x, y, area)
    @test glb.gridcontent(glb.gridcontent(ap.axis).parent).span.rows == 2:2
    @test count(b -> b isa Colorbar, gridfig.content) == 1
end
