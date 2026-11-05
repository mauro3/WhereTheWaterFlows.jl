# Catchment calculations

"""
    catchment_sinks(ctch_polygons, routing_mask, x, y[, sinks];
                    check=true, threshold_dist=Inf, connect=true)

Assign boundary sink cells to catchment polygons and reconnect fragmented labels.

Arguments
- `ctch_polygons`: vector of non-overlapping polygons, each a vector of `Point2`
  vertices or coordinate pairs.
- `routing_mask`: one-based Boolean matrix, with routed cells set to `true`.
- `x`, `y`: coordinate vectors for the matrix's first and second dimensions, in the
  same coordinate system and units as the polygons.
- `sinks`: vector of `CartesianIndex{2}` raster indices; defaults to
  `findall(get_boundary_cells(routing_mask, true))`.

KW-args
- `check=true` run a consitsency check
- `connect=true` connect disconnected parts of the different boundary groups
- `threshold_dist=Inf` distance over which cells are searched for from the polygon

Return
- `Vector{Vector{CartesianIndex{2}}}`: entry `i` contains the sink indices for
  polygon `i`, ready for `waterflows_subglacial(...; ctch_sinks=...)`.

Assign each sink to the polygon with smallest signed distance, provided that distance
is at most `threshold_dist`. By default, apply [`connect_catchment_sinks!`](@ref) as
the label cleanup. It preserves assigned sink cells and each nonempty group's
largest component; it does not smooth already-connected polygon boundaries.
Set `connect=false` to return the raw polygon assignments, including disconnected
pieces.

With `check=true`, run [`check_catchment_sinks`](@ref) to check assignment coverage.
An infinite distance threshold also requires every polygon to receive a sink before
cleanup; a finite threshold permits empty groups and unassigned distant sinks.
"""
catchment_sinks(ctch_polygons::AbstractVector, routing_mask, x, y, sinks=findall(get_boundary_cells(routing_mask, true)); kws...) =
    catchment_sinks([Point2.(cp) for cp in ctch_polygons], routing_mask, x, y, sinks; kws...)
function catchment_sinks(ctch_polygons::Vector{Vector{T}} where T<:Point2, routing_mask, x, y,
                         sinks=findall(get_boundary_cells(routing_mask, true));
                         check=true, threshold_dist=Inf, connect=true)
    ## Get all cells on the boundary of the routing_mask
    # (make two vectors: closest catchment, distance)
    closest_catch = similar(sinks, Int32)
    dist = similar(sinks, Float32)

    ## Calc signed_distance for all grounding-line pixels for all catchments
    p = Progress(length(sinks))
    Threads.@threads for ind in eachindex(sinks)
        next!(p)
        ij = sinks[ind]
        tmp = fill!(ones(length(ctch_polygons)), Inf)
        pt = Point2(x[ij[1]], y[ij[2]])

        for (i,c) in enumerate(ctch_polygons)
            tmp[i] = signed_distance(pt, c)
        end
        d, i = findmin(tmp)
        dist[ind], closest_catch[ind] = d<=threshold_dist ? (d,i) : (Inf,0)
    end
    finish!(p)
    ctch_sinks_map = fill!(similar(routing_mask, Float32), 0)
    ctch_sinks_map[sinks] .= closest_catch

    if threshold_dist==Inf
        @assert length(unique(closest_catch)) == length(ctch_polygons) "Not all ctch_polygons have boundary cells"
    end

    ## turn into ctch_sinks which can be fed to waterflows_subglacial
    ctch_sinks = [findall(ctch_sinks_map.==i) for i = 1:length(ctch_polygons)]
    if connect && connect_catchment_sinks!(ctch_sinks) > 0
        for (i, cells) in enumerate(ctch_sinks)
            ctch_sinks_map[cells] .= i
        end
    end
    check && check_catchment_sinks(ctch_sinks_map, ctch_polygons, sinks)

    return ctch_sinks #, ctch_sinks_map
end

"""
    check_catchment_sinks(ctch_sinks_map, ctch_polygons, sinks)
"""
function check_catchment_sinks(ctch_sinks_map, ctch_polygons, sinks)
    if length(unique(ctch_sinks_map)) != length(ctch_polygons)+1
        @warn "Not all ctch_polygons have sink cells"
    end
    if sum(ctch_sinks_map.>0) != length(sinks)
        @warn "Not all sinks have been assigned to catchments"
    end

    if !all(ctch_sinks_map[sinks].>0)
        @warn "Some sink cells are not part of the ctch_sinks"
    end
    if sum(ctch_sinks_map[sinks])!=sum(ctch_sinks_map)
        @warn "ctch_sinks contains more cells than there are sinks"
    end
    return nothing
end

"""
    catchment_sink_components(cells)

Return the eight-neighbour-connected components of a sink-cell group, largest first.

Arguments
- `cells` as a vector of `CartesianIndex{2}` sink-cell indices on a one-based
          raster grid, not projected coordinates.

Return
- `Vector{Vector{CartesianIndex{2}}}`, with one cell-index vector per component.

Preserve the input. Empty groups have no components; ties are broken by cell indices.
This checks connectivity, not along-boundary ordering or whether a group is unbranched.
"""
function catchment_sink_components(cells)
    remaining = Set(cells)
    components = Vector{Vector{eltype(cells)}}()
    isempty(remaining) && return components
    Iend = reduce(max, remaining)
    while !isempty(remaining)
        seed = first(remaining) # CartesianIndex minimum can be a point outside the set
        delete!(remaining, seed)
        stack = [seed]
        component = eltype(cells)[]
        while !isempty(stack)
            I = pop!(stack)
            push!(component, I)
            for J in iterate_D9(I, Iend)
                if J in remaining
                    delete!(remaining, J)
                    push!(stack, J)
                end
            end
        end
        push!(components, sort!(component; by=Tuple))
    end
    sort!(components; by=c -> (-length(c), Tuple(first(c))))
    return components
end

"""
    connect_catchment_sinks!(groups) -> Int

Reconnect fragmented sink groups in place and return the number of reassigned cells.

Arguments
- `groups`: `Vector{Vector{CartesianIndex{2}}}`, as returned by
  [`catchment_sinks`](@ref): `groups[g]` contains catchment `g`'s sink-cell indices on a
  one-based raster grid. Cell lists must be disjoint; empty catchment slots are allowed.
  These are raster indices, not projected coordinates or a raster of catchment labels.
  Mutate the inner cell vectors, retaining the outer vector and its group slots.

Return
- `Int` counts cells whose catchment index changed.

Keep the largest eight-connected component of each nonempty group fixed. Process
smaller fragments first, attaching whole fragments to an adjacent connected group.
Prefer the original label whenever it has become reachable; otherwise choose the
label with most eight-neighbour contacts, breaking ties by group index. Do not add or
remove sink cells, populate empty groups, or change an already connected assignment.
This is topological label cleanup, not independent validation of basin boundaries.

Require disjoint cell lists. Throw `ArgumentError` without modifying the input if a
fragment cannot reach any group's fixed component, for example on an isolated island
with the same label as a mainland group. Such deliberately disconnected assignments
can be retained with `catchment_sinks(...; connect=false)`.
"""
function connect_catchment_sinks!(groups)
    isempty(groups) && return 0
    labels = Dict{CartesianIndex{2},Int}(I => g for (g, cells) in enumerate(groups) for I in cells)
    length(labels) == sum(length, groups) ||
        throw(ArgumentError("Catchment sink groups must be disjoint"))
    isempty(labels) && return 0
    Iend = reduce(max, keys(labels))
    connected = Set{CartesianIndex{2}}()
    fragments = Tuple{Int,Vector{CartesianIndex{2}}}[]
    for (g, cells) in enumerate(groups)
        components = catchment_sink_components(cells)
        isempty(components) && continue
        union!(connected, first(components))
        append!(fragments, [(g, c) for c in components[2:end]])
    end
    sort!(fragments; by=f -> (length(f[2]), f[1], Tuple(first(f[2]))))
    pending = trues(length(fragments))
    changed = 0
    while any(pending)
        progress = false
        for k in eachindex(fragments)
            pending[k] || continue
            original, cells = fragments[k]
            contacts = Dict{Int,Int}()
            for I in cells, J in iterate_D9(I, Iend)
                J in connected || continue
                label = labels[J]
                contacts[label] = get(contacts, label, 0) + 1
            end
            isempty(contacts) && continue
            destination = haskey(contacts, original) ? original :
                first(sort!(collect(keys(contacts)); by=g -> (-contacts[g], g)))
            for I in cells
                labels[I] = destination
            end
            changed += destination == original ? 0 : length(cells)
            union!(connected, cells)
            pending[k] = false
            progress = true
        end
        progress || throw(ArgumentError(
            "Cannot reconnect sink groups without adding cells or moving fixed components; " *
            "use catchment_sinks(...; connect=false) for intentionally disconnected groups"))
    end
    if changed > 0
        foreach(empty!, groups)
        # Retain the raster scan order used by findall in catchment_sinks.
        for I in sort!(collect(keys(labels)); by=I -> (I[2], I[1]))
            push!(groups[labels[I]], I)
        end
    end
    return changed
end
