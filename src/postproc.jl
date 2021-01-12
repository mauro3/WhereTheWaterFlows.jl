# misc post-processing functions
"""
    prune_catchments(catchments, minsize; val=0)

Sets all catchments which are smaller than minsize to `val` (=0).
Catchments with number <1 are ignored.

Note that, as the catchments are re-numbered, the number will not correspond
to the pits anymore.
"""
function prune_catchments(catchments, minsize; val=0)
    c = copy(catchments)
    n1 = 1
    n2 = maximum(catchments)
    colormap = collect(n1:n2)
    for i=n1:n2
        if sum(x->x==i, catchments) < minsize
            colormap[i] = val
        end
    end

    # make new colors consecutive
    newcolors = unique(colormap)
    d = Dict((c=>i for (i,c) in enumerate(newcolors)))
    colormap2 = [d[c] for c in colormap]

    for i=eachindex(c)
        if c[i]>0
            c[i] = colormap2[c[i]]
        end
    end

    return c
end
