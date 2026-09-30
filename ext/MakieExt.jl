module MakieExt
# TODO:
# - the observable treatment is not done properly, i.e. some plots will not update on data updates.

# Some extension-specific threads for my information:
# - https://www.youtube.com/live/vG6ZLhe9Hns?si=7bPLZWqC-EB3AFzr&t=24578
# - https://discourse.julialang.org/t/are-extension-packages-importable/92527/18
# - https://discourse.julialang.org/t/package-extensions-structs-inside-extensions/101849
# - https://github.com/JuliaLang/Pkg.jl/pull/3552


using WhereTheWaterFlows
for fn in [:plt_dir, :plt_catchments, :plt_bnds, :plt_it, :plt_area, :plt_sinks]
    @eval import WhereTheWaterFlows:$fn
    fn!=:plt_it && @eval import WhereTheWaterFlows:$(Symbol(fn,:!))
end
const WWF=WhereTheWaterFlows
using Makie

function _tight_axis_margins!(ax=Makie.current_axis())
    if ax !== nothing
        ax.xautolimitmargin = (0.01, 0.01)
        ax.yautolimitmargin = (0.01, 0.01)
        # tightlimits!(ax, Bottom())
    end
    return nothing
end

## Preprocessing functions
sinks2inds(sinks) = ([p.I[1] for p in sinks],
                   [p.I[2] for p in sinks])
sinks2vecs(x, y, sinks) = (x[[p.I[1] for p in sinks if p!=CartesianIndex(-1,-1)]],
                           y[[p.I[2] for p in sinks if p!=CartesianIndex(-1,-1)]])


## Plotting functions

"""
    plt_dir( x, y, dir)
    
Plot `dir` as flow field. Arrow length defaults to 70% of the smaller grid
spacing (for cardinal directions); override with `lengthscale` in data units.
"""
@recipe(Plt_Dir, x, y, dir) do scene
   return Attributes(sinks=CartesianIndex{2}[], lengthscale=nothing)
end
function Makie.plot!(plot::Plt_Dir)
    (;dir, x, y, sinks) = plot
    vecfield  = lift(dir-> WWF.dir2vec.(dir, true), dir)
    vecfieldx = lift(vecfield -> [v[1] for v in vecfield], vecfield)
    vecfieldy = lift(vecfield -> [v[2] for v in vecfield], vecfield)
    lengthscale = lift(x, y, plot.lengthscale) do xs, ys, scale
        scale !== nothing && return Float64(scale)
        dx = length(xs) > 1 ? minimum(abs, diff(xs)) : Inf
        dy = length(ys) > 1 ? minimum(abs, diff(ys)) : Inf
        spacing = min(dx, dy)
        return 0.7 * (isfinite(spacing) ? spacing : 1.0)
    end
    arrows2d!(plot, x, y, vecfieldx, vecfieldy; lengthscale, align=:center)
    plt_sinks!(plot, x, y, sinks)
end

"""
    plt_area(x, y, area; prefn=log10, sinks=[], threshold=Inf, colorbar=true,
             colorbar_label="log₁₀(upstream_area)", colorbar_kwargs=(;))

Plot uparea, or another variable.

Kwargs:
- pre-proc with `prefun`, typically and by default this is `log10`
- if sinks (or pits) are passed, plot as points
- threshold the area: do not plot pixels with area below threshold
- add a colorbar on the right with `colorbar=true`
- set colorbar label with `colorbar_label`
- pass additional kwargs to `Colorbar` via `colorbar_kwargs`
"""
@recipe(Plt_Area, x, y, area) do scene
    Attributes(
        prefn = log10,
        sinks = CartesianIndex{2}[],
        threshold = Inf,
        colorbar = true,
        colorbar_label = "log₁₀(Upstream area)",
        colorbar_kwargs = (;)
    )
end
function Makie.plot!(plot::Plt_Area)
    (;x, y, area, sinks, threshold, prefn) = plot
    pl = lift(a -> prefn[].(a), area)
    if threshold[]<Inf
        pl[][area[].<threshold[]] .= NaN
    end
    plt_sinks!(plot, x, y, sinks)
    heatmap!(plot, x, y, pl)
    return plot
end

# At this stage Makie has resolved the target axis for both plt_area and plt_area!.
function Makie.plot!(ax::Makie.Axis, plot::Plt_Area)
    invoke(Makie.plot!, Tuple{Makie.AbstractAxis, Makie.AbstractPlot}, ax, plot)
    ax.xlabel = "x"
    ax.ylabel = "y"
    _tight_axis_margins!(ax)
    if plot.colorbar[]
        gc = Makie.GridLayoutBase.gridcontent(ax)
        if gc !== nothing && gc.parent !== nothing
            layout, span, side = gc.parent, gc.span, gc.side
            panel = GridLayout()
            panel[1, 1] = ax
            layout[span.rows, span.cols, side] = panel
            cbkw = plot.colorbar_kwargs[]
            hm = only(filter(p -> p isa Makie.Heatmap, plot.plots))
            Colorbar(panel[1, 2], hm;
                     label=plot.colorbar_label[], cbkw...)
        end
    end
    return plot
end

"""
    plt_sinks(x, y, sink_pits)

Plot sinks or another vector of CartesianIndices
"""
@recipe(Plt_Sinks, x, y, sinks) do scene
    Attributes()
end
function Makie.plot!(plot::Plt_Sinks)
    pp = lift((x,y,sinks) -> sinks2vecs(x, y, sinks), plot.x, plot.y, plot.sinks)
    points = lift(p -> Point2f.(p[1], p[2]), pp)
    scatter!(plot, points, color=:red, markersize=12)
end

"""
    plt_bnds(x, y, bnds)

Plot boundary points.
"""
@recipe(Plt_Bnds, x, y, bnds) do scene
    Attributes()
end
function Makie.plot!(plot::Plt_Bnds)
    (;x, y, bnds) = plot
    pp = lift((x,y,b) -> sinks2vecs(x, y, b), x, y, bnds)
    points = lift(p -> Point2f.(p[1], p[2]), pp)
    scatter!(plot, points, color=:green)
end

"""
    plt_catchments(x, y, c; minsize=0)

Plot catchments.  Catchments below `minsize` size are not plotted. With
the default `minsize=0` all catchments are plotted.

Filtering counts catchment sizes in one pass over the grid.
""" 
@recipe(Plt_Catchments, x, y, c) do scene
    Attributes(
        minsize = 0,
        colormap = :flag
        )
end
function Makie.plot!(plot::Plt_Catchments)
    (;x, y, c, minsize, colormap) = plot
    c = lift(c, minsize) do labels, threshold
        threshold > 0 ? WWF.prune_catchments(labels, threshold) : copy(labels)
    end
    colorrange = lift(labels -> (1, max(2, maximum(labels; init=0))), c)
    heatmap!(plot, x, y, c; colorrange, lowclip=(:red, 0), colormap)
    ax = Makie.current_axis()
    if ax !== nothing
        ax.xlabel = "x"
        ax.ylabel = "y"
    end
    _tight_axis_margins!()
end

# Note, this cannot be a recipe as it has several subplots
"""
    plt_it(x, y, waterflows_output, dem)

Plot DEM, uparea, flow-dir
"""
function plt_it(x, y, out::NamedTuple, dem)
    f = Makie.Figure()
    ax1 = Axis(f[1, 1]; aspect=1, title="", xautolimitmargin=(0.0, 0.0), yautolimitmargin=(0.0, 0.0))
    contour!(x, y, dem)
    tightlimits!(ax1)
    ax2 = Axis(f[2, 1]; aspect=1, title="", xautolimitmargin=(0.0, 0.0), yautolimitmargin=(0.0, 0.0))
    heatmap!(x, y, log10.(out.area))
    tightlimits!(ax2)
    ax1 = Axis(f[3, 1]; aspect=1, title="", xautolimitmargin=(0.0, 0.0), yautolimitmargin=(0.0, 0.0))
    plt_dir!(x, y, out.dir)
    tightlimits!(ax1)
    return f
end

end
