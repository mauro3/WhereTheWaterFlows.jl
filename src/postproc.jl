# Misc post-processing functions

export catchment, catchments, catchment_flux, prune_catchments, fill_dem

"""
    prune_catchments(catchments, minsize)

Sets all catchments which are smaller than minsize to `0`.
Catchments with number <1 are ignored.

Retained catchments receive consecutive positive labels.
Note that, as the catchments are re-numbered, the number will not correspond
to the sinks and pits anymore.

Counting and relabelling take O(ncells + maxlabel) time, with O(maxlabel)
auxiliary storage in addition to the output array.
"""
function prune_catchments(catchments, minsize)
    c = copy(catchments)
    n2 = max(0, maximum(catchments; init=0))
    counts = zeros(Int, n2)
    for label in catchments
        if label > 0
            counts[label] += 1
        end
    end
    colormap = zeros(Int, n2)
    nextcolor = 1
    for i=eachindex(counts)
        if counts[i] > 0 && counts[i] >= minsize
            colormap[i] = nextcolor
            nextcolor += 1
        end
    end

    for i=eachindex(c)
        if c[i]>0
            c[i] = colormap[c[i]]
        end
    end

    return c
end


"""
    catchment(dir, ij)

Calculates the catchment of one or several grid-point(s) ij.  If desired, the catchment
boundary can be calculated with `make_boundaries([c], [1])`.

Input
- dir -- direction field
- ij -- index of the point (2-Tuple, or CartesianIndex) or a Vector{CartesianIndex} for several

Returns
- catchment -- BitArray

Tip: only being off by one grid-point can make the difference
     between a tiny and a huge catchment!

See also: `catchments`
"""
catchment(dir, ij::Tuple) = catchment(dir, CartesianIndex(ij...))
function catchment(dir, ij::CartesianIndex)
    # c = fill!(similar(dir, Bool), false) # makes Matrix{Bool}
    c = fill!(similar(BitArray, axes(dir)), false) # makes BitMatrix
    # Traverse the drainage tree in up-flow direction, starting at ij.
    _catchment!(c, dir, ij)
    return c
end
function catchment(dir, ijs::Union{<:Array{CartesianIndex{2}}, CartesianIndices{2}})
    #c = fill!(similar(dir, Bool), false) # makes Matrix{Bool}
    c = fill!(similar(BitArray, axes(dir)), false) # makes BitMatrix
    for ij in ijs
        _catchment!(c, dir, ij)
    end
    return c
end
function _catchment!(c, dir, ij)
    c[ij] = true
    stack = CartesianIndex{2}[ij]
    while !isempty(stack)
        current = pop!(stack)
        for upstream in iterate_D9(current, c)
            upstream==current && continue
            c[upstream] && continue
            if flowsinto(upstream, dir[upstream], current)
                c[upstream] = true
                push!(stack, upstream)
            end
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
function catchments(dir, sinks::Union{Vector{Vector{CartesianIndex{2}}}, Vector{<:CartesianIndices{2}}};
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
        c = catchment(dir, s)
        out .+= c*i
    end
    return out
end

"""
    fill_dem(dem, sinks, dir; small=0)

Fill the pits of a DEM (apply this after applying "drainpits",
which is done by default in `waterflows`). Returns the filled DEM.

Notes:
- this is *not* needed as pre-processing step to use the flow-routing
  function `waterflows`.
- routing on a filled DEM will not produce exactly the same flow pattern: on
  the shores of lakes streams which entered the lake can now go the other way.
  I suspect on most DEMs the differences will be very minimal.
- This uses an iterative tree traversal to fill the DEM.
"""
function fill_dem(dem, sinks, dir; small=0)
    dem = copy(dem)
    Threads.@threads for sink in sinks
        _fill_ij!(-Inf, dem, sink, dir, small)
    end
    return dem
end

# The traversal goes up the catchment using `dir`. If it ever encounters a point
# which has lower elevation than the previous one, it will set that point's elevation
# and all upstream points' elevation the elevation of the first point.
function _fill_ij!(ele, dem, ij, dir, small)
    T = Union{typeof(ele),eltype(dem)}
    stack = Tuple{CartesianIndex{2},T}[(ij, ele)]
    while !isempty(stack)
        current, downstream_ele = pop!(stack)
        if downstream_ele >= dem[current]
            downstream_ele += 2*eps(downstream_ele)
            dem[current] = downstream_ele
        else
            downstream_ele = dem[current]
        end
        for upstream in iterate_D9(current, dem)
            upstream==current && continue
            if flowsinto(upstream, dir[upstream], current)
                push!(stack, (upstream, downstream_ele))
            end
        end
    end
    return nothing
end
