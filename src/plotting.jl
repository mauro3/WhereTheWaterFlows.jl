using .PyPlot
export plotarea, plotarea_dem

"""
    plotit(x, y, dem)
    plotit(x, y, waterflows_output, dem)

Plot DEM, uparea, flow-dir
"""
plotit(x, y, dem) = plotit(x, y, waterflows(dem), dem)
function plotit(x, y, waterflows_output::Tuple, dem)
    area, slen, dir, nout, nin, pits  = waterflows_output

    fig, axs = subplots(3,1)
    sca(axs[1])
    mycontour(x, y, dem)
    colorbar()
    sca(axs[2])
    heatmap(x, y, log10.(area))
    sca(axs[3])
    plotdir(x, y, dir)
end

function plotdir(x, y, dir)
    vecfield = dir2vec.(dir)
    vecfieldx = [v[1] for v in vecfield]
    vecfieldy = [v[2] for v in vecfield]
    quiver(repeat(x,1, length(y)), repeat(y,length(x),1), vecfieldx, vecfieldy)
end

pits2inds(pits) = ([p.I[1] for p in pits],
                   [p.I[2] for p in pits])
pits2vecs(x, y, pits) = (x[[p.I[1] for p in pits if p!=CartesianIndex(-1,-1)]],
                           y[[p.I[2] for p in pits if p!=CartesianIndex(-1,-1)]])

"""
    plotarea(x, y, area, pits; prefn=log10)

Plot uparea, or another variable.  If no pits should be plotted, use `[]`.
"""
function plotarea(x, y, area, pits; prefn=log10, cbar=true)
    px, py = pits2vecs(x, y, pits)
    fig, axs = subplots()
    heatmap(x, y, prefn.(area), cbar=cbar)
    scatter(px, py, 1, "r")
    fig
end

function plotarea_dem(x, y, dem, area, pits; levels=50, threshold=1/100, prefn=log10)
    fig, axs = subplots()
    mycontour(x, y, dem, levels=levels)

    threshold *= prod(size(area))
    area = copy(area)
    area[area.<threshold] .= NaN
    px, py = pits2vecs(x, y, pits)

    heatmap(x, y, prefn.(area))
    scatter(px, py, 1, "r")
end

function plotbnds(x,y,bnds)
    for b in bnds
        px, py = pits2vecs(x, y, b)
        scatter(px, py, 1, "g")
    end
end

function plotlakedepth(x, y, dem)
    area, slen, dir, nout, nin, pits  = waterflows(dem)
    demf = fill_dem(dem, pits, dir)
    heatmap(x, y, demf.-dem)
end

function heatmap(x, y, mat; cbar=true)
    dx2 = (x[2]-x[1])/2
    imshow(Array(mat'), origin="lower", extent=(x[1]-dx2,x[end]+dx2,y[1]-dx2,y[end]+dx2))
    cbar && colorbar()
end

function mycontour(x, y, mat; levels=50, cbar=true)
    contour(x, y, Array(mat'), levels, linewidths=1)
    cbar && colorbar()
end

function plot_catchments(x, y, c, lim=10)
    c = copy(c)
    tokeep = []
    for cc =1:maximum(c)
        inds = c.==cc
        if sum(inds)<lim
            c[inds] .= 0
        else
            push!(tokeep, inds)
        end
    end
    for (i,inds) in enumerate(tokeep)
        c[inds] .= i
    end
    heatmap(x,y,c)
end
