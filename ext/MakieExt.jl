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

function _tight_axis_margins!()
    ax = Makie.current_axis()
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
    
Plot `dir` as flow field
"""
@recipe(Plt_Dir, x, y, dir) do scene
   return Attributes(sinks=CartesianIndex{2}[])
end
function Makie.plot!(plot::Plt_Dir)
    (;dir, x, y, sinks) = plot
    vecfield  = lift(dir-> WWF.dir2vec.(dir, true), dir)
    vecfieldx = lift(vecfield -> [v[1] for v in vecfield], vecfield)
    vecfieldy = lift(vecfield -> [v[2] for v in vecfield], vecfield)
    arrows2d!(plot, x, y, vecfieldx, vecfieldy, lengthscale=0.005, align=:center)
    plt_sinks!(x, y, sinks)
end

"""
    plt_area(x, y, area; prefn=log10, sinks=[], threshold=Inf, colorbar=true,
             colorbar_label="log10(upstream_area)", colorbar_kwargs=(;))

Plot uparea, or another variable.

Kwargs:
- pre-proc with `prefun`, typically and by default this is `log10`
- if sinks (or pits) are passed, plot as points
- threshold the area: do not plot pixels with area below threshold
- add a colorbar with `colorbar=true`
- set colorbar label with `colorbar_label`
- pass additional kwargs to `Colorbar` via `colorbar_kwargs`
"""
@recipe(Plt_Area, x, y, area) do scene
    Attributes(
        prefn = log10,
        sinks = CartesianIndex{2}[],
        threshold = Inf,
        colorbar = true,
        colorbar_label = "log10(Upstream area)",
        colorbar_kwargs = (;)
    )
end
function Makie.plot!(plot::Plt_Area)
    (;x, y, area, sinks, threshold, prefn, colorbar, colorbar_label, colorbar_kwargs) = plot
    pl = lift(a -> prefn[].(a), area)
    if threshold[]<Inf
        pl[][area[].<threshold[]] .= NaN
    end
    plt_sinks!(plot, x, y, sinks)
    hm = heatmap!(plot, x, y, pl)
    ax = Makie.current_axis()
    if ax !== nothing
        ax.xlabel = "x"
        ax.ylabel = "y"
    end
    if colorbar[]
        fig = Makie.current_figure()
        if fig !== nothing
            cbkw = colorbar_kwargs[]
            if haskey(cbkw, :label)
                Colorbar(fig[:, end+1], hm; cbkw...)
            else
                Colorbar(fig[:, end+1], hm; label=colorbar_label[], cbkw...)
            end
        end
    end
    _tight_axis_margins!()
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
    px, py = lift(x->x, pp[])
    scatter!(plot, px, py, color=:red, markersize=12)
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
    px, py = lift((x,y,b) -> sinks2vecs(x, y, b), x, y, bnds)
    scatter!(px, py, color=:green)
end

"""
    plt_catchments(x, y, c; minsize=0)

Plot catchments.  Catchments below `minsize` size are not plotted. With
the default `minsize=0` all catchments are plotted.

Note, `minsize>0` can be quite slow to compute.
""" 
@recipe(Plt_Catchments, x, y, c) do scene
    Attributes(
        minsize = 0,
        colormap = :flag
        )
end
function Makie.plot!(plot::Plt_Catchments)
    (;x, y, c, minsize, colormap) = plot
    c = lift(copy, c)
    tokeep = BitMatrix[]
    if minsize[]>0
        for cc = 1:maximum(c[])
            inds = c[].==cc
            if sum(inds)<minsize[]
                c[][inds] .= 0
            else
                push!(tokeep, inds)
            end
        end
        for (i,inds) in enumerate(tokeep)
            c[][inds] .= i
        end
    end
    heatmap!(plot, x, y, c; colorrange=(1,maximum(c[])), lowclip=(:red, 0), colormap)
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
