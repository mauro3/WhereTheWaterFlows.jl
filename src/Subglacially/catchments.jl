# Catchment calculations

const CI = CartesianIndex{2}

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
largest D4 component, preferring strong connections with diagonal fallback.
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
    if connect
        n = connect_catchment_sinks!(ctch_sinks)
        if check && n > 0
            for (i, cells) in enumerate(ctch_sinks)
                ctch_sinks_map[cells] .= i
            end
        end
    end
    check && check_catchment_sinks(ctch_sinks_map, ctch_polygons, sinks)

    return ctch_sinks
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
    catchment_sink_components(cells; diagonal=true)

Return connected components of a sink-cell group, largest first.

Use eight-neighbour connectivity by default, or D4 with `diagonal=false`.

Arguments
- `cells` as a vector of `CartesianIndex{2}` sink-cell indices

Return
- `Vector{Vector{CartesianIndex{2}}}`, with one cell-index vector per component.

Preserve the input. Empty groups have no components; ties are broken by cell indices.
This checks connectivity, not along-boundary ordering or whether a group is unbranched.
"""
function catchment_sink_components(cells; diagonal=true)
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
                !diagonal && sum(abs, Tuple(J-I)) != 1 && continue
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

Keep the largest D4 component of each nonempty group fixed. Process smaller
fragments first, exhausting D4 attachments before allowing a diagonal attachment.
After each diagonal attachment, retry D4 attachments. Within each connectivity tier,
prefer the original label, then most contacts, then lowest group index. D8-only
appendices remain possible where the underlying boundary has no D4 connection.
Do not add or remove sink cells or populate empty groups. Cells in an originally
D4-connected group keep their label. This is topological label cleanup, not
independent validation of basin boundaries.

Require disjoint cell lists. Throw `ArgumentError` without modifying the input if a
fragment cannot reach any group's fixed component, for example on an isolated island
with the same label as a mainland group. Such deliberately disconnected assignments
can be retained with `catchment_sinks(...; connect=false)`.
"""
function connect_catchment_sinks!(groups)
    isempty(groups) && return 0
    # label each cell's group. Of form: CartesianIndex(12, 34) => 5
    labels = Dict{CI,Int}(I => g for (g, cells) in enumerate(groups) for I in cells)
    length(labels) != sum(length, groups) &&
        throw(ArgumentError("Catchment sink groups must be disjoint"))
    isempty(labels) && return 0
    Iend = reduce(max, keys(labels))
    # store cells connected to main bit
    connected = Set{CI}()
    # and disconnect bits
    fragments = Tuple{Int,Vector{CI}}[]
    # 1. Keep one fixed piece per group
    for (g, cells) in enumerate(groups)
        components = catchment_sink_components(cells; diagonal=false)
        isempty(components) && continue
        union!(connected, first(components))
        append!(fragments, [(g, c) for c in components[2:end]])
    end
    # 2. Exhaust strong attachments before taking one weak step and trying again.
    sort!(fragments; by=f -> (length(f[2]), f[1], Tuple(first(f[2]))))
    pending = trues(length(fragments))
    changed = 0
    while any(pending)
        progress = false
        for diagonal in (false, true)
            for k in eachindex(fragments)
                !pending[k] && continue
                original, cells = fragments[k]
                contacts = Dict{Int,Int}()
                for I in cells, J in iterate_D9(I, Iend)
                    !diagonal && sum(abs, Tuple(J-I)) != 1 && continue
                    !(J in connected) && continue
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
                diagonal && break
            end
            progress && break
        end
        !progress && throw(ArgumentError(
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

# Circular order in (dimension 1, dimension 2), not raster-scan order.
const D8_POS = (
    CI(1, 0), CI(1, 1), CI(0, 1),
    CI(-1, 1), CI(-1, 0), CI(-1, -1),
    CI(0, -1), CI(1, -1),
)
const D4_POS = D8_POS[1:2:8]

# Determine if I is in region :inside (within mask minus gl-cells) or
# in :outside (where mask==false). I outside size(mask) returns false.
function in_region(I::CI, g::Set, mask, region::Symbol)
    checkbounds(Bool, mask, I) || return false
    return region == :inside ? mask[I] && !(I in g) : !mask[I]
end

"""
    touches_side(I, J, g, mask, region, side) -> Bool

Test whether `region` touches `side` of the directed step `I -> J`.

Supply distinct D8-neighbouring 2D indices `I` and `J`, a Boolean routing `mask`,
and the complete GL cell set `g`, including other outlets' GL cells. Use
`region=:inside` for mask==true excluding `g`, or `:outside` for mask==false cells.
Out-of-bounds cells satisfy neither region. Choose `side=:left` or `:right`.

WWF treats the first matrix dimension (rows) as x and the second (columns) as y.
Left is defined in that Cartesian frame: for `d = J - I = (dx, dy)`, the left
normal is `(-dy, dx)`; right uses the opposite normal, not ordinary matrix-display
orientation. For example, a step `(1, 0)` has its left side towards `(0, 1)`.

For a non-diagonal step, test the cells offset by that normal from both endpoints;
return `true` if either is in `region`, not necessarily both. For a diagonal step,
test the single shared cardinal neighbour on the requested side. Only these
side-neighbours are classified; membership of `I` and `J` in `g` is assumed but not checked.

See also [`in_region`](@ref), [`step_gl`](@ref).
"""
function touches_side(I::CI, J::CI, g::Set, mask, region::Symbol, side::Symbol)
    d = J - I
    sign = side == :left ? 1 : -1
    # a non-diagonal step
    if d[1] == 0 || d[2] == 0
        normal = sign * CI(-d[2], d[1])
        return in_region(I + normal, g, mask, region) ||
               in_region(J + normal, g, mask, region)
    end
    # A diagonal step has one shared neighbour on each side.
    K = I + (sign * d[1] * d[2] > 0 ? CI(0, d[2]) :
                                              CI(d[1], 0))
    return in_region(K, g, mask, region)
end

"""
    step_gl(previous, I, g, mask; region=:inside, side=:left, visited=(), boundary=g)

Choose a neighbour relative to the arrival direction, or return `nothing`.

Inside steps are D4 and must touch interior on `side`. Outside steps are D8;
prefer the outside on `side`, then allow remaining branches for complete pickup.
Pass the allowed GL cells as a set `g`; `visited` excludes cells already recorded.
Use `boundary` for the complete GL set when walking just one outlet group.
`previous` must be a D4 neighbour of `I` for inside steps, or a D8 neighbour otherwise.
"""
function step_gl(previous::CI, I::CI, g::Set, mask; region=:inside, side=:left, visited=(), boundary=g)
    offsets = region == :inside ? D4_POS : D8_POS
    back = findfirst(==(previous - I), offsets)
    turn = side == :left ? -1 : 1
    for hugging in (true, false)
        region == :inside && !hugging && break
        for n in 1:length(offsets)
            J = I + offsets[mod1(back + turn * n, length(offsets))]
            J in g && !(J in visited) || continue
            touches_side(I, J, boundary, mask, region, side) == hugging || continue
            return J
        end
    end
    return nothing
end

"""
    trace_inside(g, mask, start, side, previous; traversed=nothing, boundary=g)

Trace a D4 GL path from `start`, keeping interior on `side`.

Use `previous` to define the arrival direction and `boundary` to classify interior.
Stop at a dead end, a U-turn, or before repeating a directed edge. Return the cell
sequence including `start`; cells may repeat. Optionally add edges to `traversed`,
oriented with interior on the left.
"""
function trace_inside(g::Set, mask, start::CI, side::Symbol, previous::CI;
                      traversed=nothing, boundary=g)
    path = [start]
    steps = Set{Tuple{CI, CI}}()
    I = start
    while true
        J = step_gl(previous, I, g, mask; side, boundary)
        (J === nothing || J == previous || (I, J) in steps) && break
        push!(steps, (I, J))
        # Store both trace directions with interior on the left.
        traversed === nothing || push!(traversed, side == :left ? (I, J) : (J, I))
        push!(path, J)
        previous, I = I, J
    end
    return path
end

"""
    walk_gl_inside(gl, mask, start=first(gl); boundary=Set(gl))

Follow the local D4 inside walk in both directions, recording each cell once.

Supply nonempty GL cell indices `gl` and a Boolean routing `mask` with active cells
set to `true`. Choose a `start` in `gl` touching interior. Interior is an active cell
outside `boundary`; supply the complete GL set when `gl` is only one outlet group.
This local walk may omit appendices and other inside walks.

See also [`walk_gl`](@ref) to discover all walks or select the longest.
"""
function walk_gl_inside(gl::AbstractVector, mask, start=first(gl); boundary=Set(gl))
    # Do
    # 
    g = Set(gl)
    start in g || error("Start is not on the grounding line: $start")
    any(d -> in_region(start + d, boundary, mask, :inside), D8_POS) ||
        error("Start does not touch interior: $start")
    # Choose a real seed edge before tracing, so the incoming cell is never invented.
    previous = start + first(D4_POS)
    side = :left
    J = step_gl(previous, start, g, mask; side, boundary)
    if J === nothing
        side = :right
        J = step_gl(previous, start, g, mask; side, boundary)
    end
    J === nothing && return [start]
    forward = trace_inside(g, mask, J, side, start; boundary)
    backward = trace_inside(g, mask, start, side == :left ? :right : :left, J; boundary)
    return unique([reverse(backward); forward])
end

"""
    walk_gl_outside(gl, mask, start=first(gl))

Visit every cell of a D8-connected GL once, following the outside first.

Supply nonempty GL cell indices `gl` and a Boolean routing `mask`, with `false`
cells marking exterior. Start at a cell in `gl`; throw if it cannot reach all cells.
Use a stack to return from branches without recording cells twice. Consecutive
output cells can therefore jump when a branch has been exhausted.

See also [`walk_gl`](@ref), [`walk_gl_inside`](@ref).
"""
function walk_gl_outside(gl, mask, start=first(gl))
    g = Set(gl)
    start in g || error("Start is not on the grounding line: $start")
    path, visited = [start], Set([start])
    stack = [(start + first(D8_POS), start)]
    while !isempty(stack)
        previous, I = last(stack)
        J = step_gl(previous, I, g, mask; region=:outside, side=:right, visited)
        if J === nothing
            pop!(stack)
        else
            push!(path, J)
            push!(visited, J)
            push!(stack, (I, J))
        end
    end
    length(path) == length(g) || error("GL is not D8-connected from $start")
    return path
end

"""
    walk_gl(gl, mask, start=first(gl); cull_inside=true, boundary=Set(gl))

Return inside walks and the complete outside ordering as `(; inside, outside)`.

Supply nonempty, D8-connected GL cell indices `gl` and a Boolean routing `mask`.
Active non-GL cells are interior; inactive cells are exterior. Supply the complete GL
set as `boundary` when `gl` contains only one outlet group, so other groups' GL cells
are not classified as interior. Neither the mask nor the cell lists are changed.

By default, `inside` is the longest D4 walk as a vector of unique cells (first on
ties). With `cull_inside=false`, return all discovered walks as separate vectors;
walks may share junctions. Return an empty `inside` if no GL cell touches interior.
Inside walks may omit appendices. Discover them along the outside ordering even
when `start` is on an appendix or a smaller inside walk.

The `outside` vector visits every GL cell once, starting at `start`. Consecutive
output cells may jump after unrecorded backtracking; it is not a continuous polyline.

See also [`walk_gl_inside`](@ref), [`walk_gl_outside`](@ref).
"""
function walk_gl(gl, mask, start=first(gl); cull_inside=true, boundary=Set(gl))
    outside = walk_gl_outside(gl, mask, start)
    g = Set(gl)
    walks = Vector{CI}[]
    traversed = Set{Tuple{CI, CI}}()
    for I in outside
        any(d -> in_region(I + d, boundary, mask, :inside), D8_POS) || continue
        for d in D4_POS
            J = I + d
            J in g && !((I, J) in traversed) || continue
            touches_side(I, J, boundary, mask, :inside, :left) || continue
            push!(traversed, (I, J))
            forward = trace_inside(g, mask, J, :left, I; traversed, boundary)
            backward = trace_inside(g, mask, I, :right, J; traversed, boundary)
            path = unique([reverse(backward); forward])
            push!(walks, path)
        end
    end
    # Retain inside-touching cells with no inside step, without duplicating endpoints.
    covered = Set(Iterators.flatten(walks))
    for I in outside
        I in covered && continue
        any(d -> in_region(I + d, boundary, mask, :inside), D8_POS) || continue
        push!(walks, [I])
    end
    inside = cull_inside ? (isempty(walks) ? CI[] : argmax(length, walks)) : walks
    return (; inside, outside)
end
