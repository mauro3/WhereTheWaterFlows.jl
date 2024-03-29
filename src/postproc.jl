# Misc post-processing functions

export catchment, catchments, catchment_flux, prune_catchments, fill_dem

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


"""
    catchment(dir, ij, dem=nothing)

Calculates the catchment of one or several grid-point ij.  If desired, its
boundary can be calculated with `make_boundaries([c], [1])`.

Input
- dir -- direction field
- ij -- index of the point (2-Tuple, or CartesianIndex) or a Vector{CartesianIndex} for several
- dem -- if supplied, then its NaN-points will not be considered. NOTE: this slow down the calculation
         by potentially a large amount!

Returns
- catchment -- BitArray

Tip: only being off by one grid-point can make the difference
     between a tiny and a huge catchment!

See also: `catchments`
"""
catchment(dir, ij::Tuple, dem=nothing) = catchment(dir, CartesianIndex(ij...), dem)
function catchment(dir, ij::CartesianIndex, dem=nothing)
    # c = fill!(similar(dir, Bool), false) # makes Matrix{Bool}
    c = fill!(similar(BitArray, axes(dir)), false) # makes BitMatrix
    # recursively traverse the drainage tree in up-flow direction,
    # starting at ij
    _catchment!(c, dir, ij)
    if dem!==nothing
        # TODO: this is very non-performant
        c[isnan.(dem)] .= false
    end
    return c
end
function catchment(dir, ijs::Union{<:Array{CartesianIndex{2}}, CartesianIndices{2}}, dem=nothing)
    #c = fill!(similar(dir, Bool), false) # makes Matrix{Bool}
    c = fill!(similar(BitArray, axes(dir)), false) # makes BitMatrix
    for ij in ijs
        _catchment!(c, dir, ij)
    end
    if dem!==nothing
        # TODO: this is very non-performant
        c[isnan.(dem)] .= false
    end
    return c
end
function _catchment!(c, dir, ij)
    c[ij] = true
    # proc upstream points
    for IJ in iterate_D9(ij, c)
        ij==IJ && continue
        c[IJ] && continue # this hits a point which is already in the catchment,
                          # thus nothing more to do
        if flowsinto(IJ, dir[IJ], ij)
            _catchment!(c, dir, IJ)
        end
    end
    return nothing
end

"""
    catchment_flux(cellarea, c, color) = sum(cellarea[c.==color])
    catchment_flux(cellarea, c::Union{BitArray, Matrix{Bool}})

The total flux, i.e. input, in one catchment.
"""
catchment_flux(cellarea, c, color) = sum(cellarea[c.==color])
catchment_flux(cellarea, c::Union{BitArray, AbstractMatrix{<:Bool}}) = sum(cellarea[c])

"""
    catchments(dir, sinks::Union{Vector{Vector{CartesianIndex{2}}}, Vector{<:CartesianIndices{2}}}, dem=nothing;
                    check_sinks_overlap=true)

Make a map of catchments from different (non-overlapping) sinks.

See also: `catchment`
"""
function catchments(dir, sinks::Union{Vector{Vector{CartesianIndex{2}}}, Vector{<:CartesianIndices{2}}},
                    dem=nothing;
                    check_sinks_overlap=true)

    ncs = length(sinks)
    if check_sinks_overlap
        for s in sinks
            for loc in s
                for ss in sinks
                    ss===s && continue
                    loc in ss && error("Detected overlapping skinks.")
                end
            end
        end
    end
    @assert ncs<255 "More than 255 sinks not supported (yet)." # then use something else than UInt8
    out = fill!(similar(dir, UInt8), 0)
    for (i,s) in enumerate(sinks)
        c = catchment(dir, s, dem)
        out += c*i
    end
    return out
end

"""
    fill_dem(dem, pits, dir)

Fill the pits (aka sinks) of a DEM (apply this after applying "drainpits",
which is the default). Returns the filled DEM.

Note, this is *not* needed as pre-processing step to use the flow-routing
function `waterflows`.

This uses a tree traversal to fill the DEM. It does it depth-first (as it
is easier) which may lead to a stack overflow on a large DEM.
"""
function fill_dem(dem, pits, dir; small=0.0)
    dem = copy(dem)
    Threads.@threads for pit in pits
        _fill_ij!(-Inf, dem, pit, dir, small)
    end
    return dem
end

# The recursion goes up the catchment using `dir`.  If it ever encounters a point
# which has lower elevation than the previous one, it will set that point's elevation
# and all upstream points' elevation the elevation of the first point.
function _fill_ij!(ele, dem, ij, dir, small)
    if ele > dem[ij]
        dem[ij] = ele
        ele += small
    else
        ele = dem[ij]
    end
    # proc upstream points
    for IJ in iterate_D9(ij, dem)
        ij==IJ && continue
        if flowsinto(IJ, dir[IJ], ij)
            _fill_ij!(ele, dem, IJ, dir, small)
        end
    end
end
