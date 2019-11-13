
using .PyPlot

"""
    plotit(xs, ys, dem)

Plot DEM, uparea, flow-dir
"""
function plotit(xs, ys, dem)
    dx2 = step(xs)/2
    area, slen, dir, nout, nin, pits  = waterflows(dem)

    fig, axs = subplots(3,1)
    sca(axs[1])
    contour(xs, ys, Array(dem'))
    colorbar()
    sca(axs[2])
    heatmap(xs, ys, log10.(area))
    sca(axs[3])
    plotdir(xs, ys, dir)
end

function plotdir(xs, ys, dir)
    vecfield = dir2vec_0.(dir)
    vecfieldx = [v[1] for v in vecfield]
    vecfieldy = [v[2] for v in vecfield]
    quiver(repeat(xs,1, length(ys)), repeat(ys,length(xs),1), vecfieldx, vecfieldy)
end

pits2inds(pits) = ([p.I[1] for p in pits],
                   [p.I[2] for p in pits])
pits2vecs(xs, ys, pits) = (xs[[p.I[1] for p in pits if p!=CartesianIndex(-1,-1)]],
                           ys[[p.I[2] for p in pits if p!=CartesianIndex(-1,-1)]])

"""
    plotarea(xs, ys, dem; prefn=log10)
    plotarea(xs, ys, area, pits; prefn=log10)

Plot uparea
"""
function plotarea(xs, ys, dem; prefn=log10)
    area, slen, dir, nout, nin, pits  = waterflows(dem)
    plotarea(xs, ys, area, pits, prefn=prefn)
end
function plotarea(xs, ys, area, pits; prefn=log10)
    px, py = pits2vecs(xs, ys, pits)
    fig, axs = subplots()
    heatmap(xs, ys, prefn.(area))
    scatter(px, py, 1, "r")
end

function plotlakedepth(x, y, dem)
    area, slen, dir, nout, nin, pits  = waterflows(dem)
    demf = fill_dem(dem, pits, dir)
    heatmap(x, y, demf.-dem)
end

function heatmap(xs, ys, mat)
    dx2 = step(xs)/2
    imshow(Array(mat'), origin="lower", extent=(xs[1]-dx2,xs[end]+dx2,ys[1]-dx2,ys[end]+dx2))
    colorbar()
end
