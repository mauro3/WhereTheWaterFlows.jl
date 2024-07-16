module MakieExt
# TODO:
# - the observable treatment is not done properly, i.e. some plots will not update on data updates.

# Some extension-specific threads for my information:
# - https://www.youtube.com/live/vG6ZLhe9Hns?si=7bPLZWqC-EB3AFzr&t=24578
# - https://discourse.julialang.org/t/are-extension-packages-importable/92527/18
# - https://discourse.julialang.org/t/package-extensions-structs-inside-extensions/101849
# - https://github.com/JuliaLang/Pkg.jl/pull/3552


using WhereTheWaterFlows
for fn in [:plt_dir, :plt_catchments, :plt_lakedepth, :plt_bnds, :plt_it, :plt_area, :plt_sinks] 
    @eval import WhereTheWaterFlows:$fn
end
const WWF=WhereTheWaterFlows
using Makie

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
    arrows!(plot, x, y, vecfieldx, vecfieldy, lengthscale=0.005, align=:center)
    plt_sinks!(x, y, sinks)
end

"""
    plt_area(x, y, area; prefn=log10, sinks=[], threshold=Inf)

Plot uparea, or another variable.

Kwargs:
- pre-proc with `prefun`, typically this is `log10`
- if sinks (or pits) are passed, plot as points
- threshold the area
"""
@recipe(Plt_Area, x, y, area) do scene
    Attributes(
        prefn = log10,
        sinks = CartesianIndex{2}[],
        threshold = Inf
    )
end
function Makie.plot!(plot::Plt_Area)
    (;x, y, area, sinks, threshold) = plot
    prefn = plot.prefn[]
    pl = lift(a -> prefn.(a), area)
    if threshold[]<Inf
        pl[][area[].<threshold[]] .= NaN
    end
    plt_sinks!(plot, x, y, sinks)
    heatmap!(plot, x, y, pl)
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
    scatter!(plot, px, py, color=:red, markersize=25)
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
    px, py = lift((x,y,b) -> sink_pits2vecs(x, y, b), x, y, bnds)
    scatter!(px, py, color=:green)
end

"""
    plt_catchments(x, y, c; minsize=10)

Plot catchments.  Catchments below `minsize` size are not plotted.
""" 
@recipe(Plt_Catchments, x, y, c) do scene
    Attributes(
        minsize = 10,
        colormap = :flag
        )
end
function Makie.plot!(plot::Plt_Catchments)
    (;x, y, c, minsize, colormap) = plot
    c = lift(copy, c)
    tokeep = BitMatrix[]
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
    heatmap!(plot, x, y, c; colorrange=(1,maximum(c[])), lowclip=(:red, 0), colormap)
end

"""
    plt_lakedepth(x, y, dem; lowerlimit=0)
    plt_lakedepth(x, y, dem, dir, sinks; lowerlimit=0)

Plot lake depth.  
"""
@recipe(Plt_Lakedepth, x, y, dem, dir, sinks) do scene
    Attributes(
        lowerlimit=0
    )
end
function Makie.plot!(plot::Plt_Lakedepth)
    (;x, y, dem, dir, sinks, lowerlimit) = plot
    lakedepth = WWF.fill_dem(dem[], sinks[], dir[]) .- dem[]
    heatmap!(plot, x, y, lakedepth; colorrange=(lowerlimit[], maximum(lakedepth)), lowclip=(:red, 0))
    # TODO: colorbar
end
# specialize for 3-arg call:
const _Plt_Lakedepth_type = Plt_Lakedepth{<:NTuple{3,Any}}
argument_names(::Type{_Plt_Lakedepth_type}) = (:x, :y, :dem)
function Makie.plot!(plot::_Plt_Lakedepth_type)
    (;x, y, dem, lowerlimit) = plot
    (;lowerlimit) = plot
    dir, sinks   = waterflows(dem[])[[3,6]]
    plt_lakedepth!(plot, x, y, dem, dir, sinks; lowerlimit)
end

# Note, this cannot be a recipe as it has several subplots
"""
    plt_it(x, y, dem)
    plt_it(x, y, waterflows_output, dem)

Plot DEM, uparea, flow-dir
"""
plt_it(x, y, dem) = plt_it(x, y, waterflows(dem), dem)
function plt_it(x, y, waterflows_output::Tuple, dem)
    area, slen, dir, nout, nin, sinks, sink_pits  = waterflows_output

    f = Makie.Figure()
    ax1 = Axis(f[1, 1]; aspect=1, title="")
    contour!(x, y, dem)
    ax2 = Axis(f[2, 1]; aspect=1, title="")
    heatmap!(x, y, log10.(area))
    ax1 = Axis(f[3, 1]; aspect=1, title="")
    plt_dir!(x, y, dir)
    return f
end

end
