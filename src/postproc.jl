# Misc post-processing functions

export catchment, catchment!, catchments, catchment_flux, prune_catchments, fill_dem
public last_visit_fluxes

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

See also: [`catchment!`](@ref), `catchments`.
"""
function catchment(dir, ij)
    out = similar(BitArray, axes(dir))
    catchment!(out, dir, ij)
    return out
end

"""
    catchment!(mask, dir, seeds) -> mask

Clear and fill a Boolean `mask` with the upstream catchment of `seeds`.

- Supply `mask` with the same axes as `dir`
- `seeds` a coordinate tuple, a `CartesianIndex`, or a collection of Cartesian indices.

Clears the mask once per call, but not between seeds; repeated seeds are ignored.

See also [`catchment`](@ref).
"""
catchment!(mask::AbstractMatrix{Bool}, dir, ij::Tuple{Integer,Integer}) =
    catchment!(mask, dir, CartesianIndex(ij...))
catchment!(mask::AbstractMatrix{Bool}, dir, ij::CartesianIndex) =
    catchment!(mask, dir, (ij,))
function catchment!(mask::AbstractMatrix{Bool}, dir, seeds)
    axes(mask) == axes(dir) || throw(DimensionMismatch("Mask and direction axes differ"))
    # Expose packed-bit storage without changing the output container.
    # https://github.com/rafaqz/DimensionalData.jl/issues/1244
    fill!(vec(mask), false)
    for ij in seeds
        _catchment!(mask, dir, ij)
    end
    return nothing
end
function _catchment!(c, dir, ij)
    c[ij] && return nothing
    c[ij] = true
    # proc upstream points
    for IJ in iterate_D9(ij, c)
        ij==IJ && continue
        c[IJ] && continue # this hits a point which is already processed,
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
catchment_flux(cellarea, c, color) = catchment_flux(cellarea, c.==color)
function catchment_flux(cellarea, c::Union{BitArray, AbstractMatrix{<:Bool}})
    size(cellarea) == size(c) || throw(DimensionMismatch("Flux and mask sizes differ"))
    total = zero(Base.promote_op(Base.add_sum, eltype(cellarea), eltype(cellarea)))
    for i in Base.LogicalIndex{Int}(vec(c))
      total += cellarea[i]
    end
    return total
end

"""
    last_visit_fluxes(dir, terminals, labels, source; nlabels=nothing)

Allocate each labelled source to its last same-label cell downstream.

Supply a draining D8 forest, its unique terminals, nonnegative integer region
`labels`, and additive source weights on the same one-based grid. Label zero is
ignored. Regions may be disconnected and paths may leave and re-enter a region.
A positive-labelled terminal acts as a last-visit location at that terminal.
All intended source cells must drain to one of the supplied terminals.

Return `cells`, their region `labels`, accumulated `flux`, total `source_flux`
per region, and `visited_cells`. Include zero-flux last-visit cells. Specify
`nlabels` to retain absent region IDs in `source_flux`. No feedback or routing
changes are applied; this allocates the supplied source on existing directions.

Traverse each tree once, using a downstream-most anchor per label on the current
path. Restore anchors when returning from branches. Storage is one grid of source
precision, a label-sized anchor table, sparse outputs and the recursion stack.
"""
function last_visit_fluxes(dir::AbstractMatrix, terminals,
                           labels::AbstractMatrix{<:Integer}, source::AbstractMatrix;
                           nlabels=nothing)
    Base.require_one_based_indexing(dir, labels, source)
    axes(dir) == axes(labels) == axes(source) ||
        throw(DimensionMismatch("Directions, labels and source grids differ"))
    lo, hi = extrema(labels)
    n = isnothing(nlabels) ? Int(hi) : Int(nlabels)
    0 <= lo <= hi <= n || throw(ArgumentError("Labels must lie between zero and nlabels"))
    T = float(eltype(source))
    storage = (;dir, labels, source, anchors=zeros(Int,n),
               flux=zeros(T,size(dir)), source_flux=zeros(T,n),
               linear=LinearIndices(dir), cells=CartesianIndex{2}[], regions=Int[])
    visited = 0
    for terminal in terminals
        visited += _last_visit_fluxes!(storage, terminal)
    end
    return (;cells=storage.cells, labels=storage.regions,
            flux=storage.flux[storage.cells], source_flux=storage.source_flux,
            visited_cells=visited)
end

function _last_visit_fluxes!(s, I)
    # The first occurrence while walking upstream is the LAST one downstream.
    # Keep it for descendants, including paths that leave and re-enter the label.
    @inbounds label = Int(s.labels[I])
    previous = 0
    if label > 0
        @inbounds previous = s.anchors[label]
        if previous == 0
            @inbounds s.anchors[label] = s.linear[I]
            push!(s.cells,I)
            push!(s.regions,label)
        end
        @inbounds begin
            value = s.source[I]
            s.flux[s.anchors[label]] += value
            s.source_flux[label] += value
        end
    end
    visited = 1
    for J in iterate_D9(I,s.dir)
        I == J && continue
        if flowsinto(J,s.dir[J],I)
            visited += _last_visit_fluxes!(s,J)
        end
    end
    # A sibling branch must never inherit an anchor from the branch just finished.
    label > 0 && (@inbounds s.anchors[label] = previous)
    return visited
end

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
- This uses a tree traversal to fill the DEM. It does it depth-first (as it
  is easier) which may lead to a stack overflow on a large DEM.
"""
function fill_dem(dem, sinks, dir; small=0)
    dem = copy(dem)
    Threads.@threads for sink in sinks
        _fill_ij!(-Inf, dem, sink, dir, small)
    end
    return dem
end

# The recursion goes up the catchment using `dir`.  If it ever encounters a point
# which has lower elevation than the previous one, it will set that point's elevation
# and all upstream points' elevation the elevation of the first point.
function _fill_ij!(ele, dem, ij, dir, small)
    if ele >= dem[ij]
        ele += 2*eps(ele)
        dem[ij] = ele
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
    return nothing
end
