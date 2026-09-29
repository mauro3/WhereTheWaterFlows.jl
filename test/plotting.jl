using Test, WhereTheWaterFlows
using CairoMakie

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
