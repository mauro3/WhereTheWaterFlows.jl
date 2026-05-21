#########
# Helpers
#########

"""
    boxcar(A::AbstractArray, window::AbstractArray, weights)

Moving average filter with spatially dependent filter window and weighting of cells.
Filters over ±window.
"""
function boxcar(A::AbstractArray{T,N}, window::AbstractArray{<:Integer,N},
                weights::AbstractArray{<:Number,N}=ones(size(A)...)) where {T,N}
    out = similar(A)
    R = CartesianIndices(size(A))
    I1, Iend = first(R), last(R)
    Threads.@threads for I in R # @inbounds does not help
        out[I] = NaN
        n, s = 0, zero(eltype(out))
        I_ul = CartesianIndex(I1.I.*window[I])
        for J in CartesianIndices(UnitRange.(max(I1, I-I_ul).I , min(Iend, I+I_ul).I) )
            # used to be CartesianRange(max(I1, I-I_l), min(Iend, I+I_u) )
            # now it is probably something simpler than what I use above
            AJ, w = A[J], weights[J]
            if !isnan(AJ) && w!=0
                s += A[J]
                n += 1
            end
        end
        if n==0
            out[I] = A[I]
        else
            out[I] = s/n # note: ==NaN if n==s==0
        end
    end
    out
end

"""
    smooth_surface(x, y, surfdem, beddem, icethicknesses, mask= surfdem.>=beddem; minwindow=0)

Smooth the surface DEM of ±`icethicknesses`.  Only smoothes where
there is ice and only using ice-covered cells.
"""
function smooth_surface(x, y, surfdem, beddem, icethicknesses, mask=surfdem.>=beddem; minwindow=0)
    dx = x[2]-x[1]
    @assert y[2]-y[1]==dx
    thickness = surfdem .- beddem
    thickness[isnan.(thickness)] .= 0
    window = round.(Int, icethicknesses .* thickness ./ dx)
    window[window.<minwindow] .= minwindow
    surface_smooth = boxcar(surfdem, 
    window, mask)
    return surface_smooth
end

"Return CartesianIndices corresponding to the 24 neighbors and the point itself."
function iterate_D25(I::CartesianIndex, ar::AbstractMatrix)
    R = CartesianIndices(ar)
    two = CartesianIndex(2,2)
    I1, Iend = first(R), last(R)
    return max(I1, I-two):min(Iend, I+two)
end

"""
    mask_contiguous(mask, IJ, out=similar(mask, Bool))

Create a new mask (or update the optional 3rd argument) which
- marks all points connected to IJ and
- where each masked point has the same value as mask[IJ]
"""
function mask_contiguous(mask, IJ, out=similar(mask, Bool))
    fill!(out, false)
    val = mask[IJ]
    queue = [IJ]
    out[IJ] = true
    while length(queue)>0
        # get first element
        IJ = popfirst!(queue)
        # push onto queue
        for ij in iterate_D9(IJ, out)
            out[ij] && continue # already done or queued
            if mask[ij]==val
                out[ij] = true
                push!(queue, ij) # queue new point
            end
        end
        # make sure it does not blow up
        length(queue)>prod(size(mask)) && error("oh no, queue got too big")
    end
    return out
end

"""
    get_boudary_cells(mask, value)

Get all the cells of mask of value `value` which border
on cells of different value.
"""
function get_boudary_cells(mask, value)
    out = similar(mask, Bool) .* false

    for IJ in CartesianIndices(mask)
        mask[IJ]!=value && continue # don't process
        for ij in iterate_D9(IJ, mask)
            ij==IJ && continue # don't process the point itself
            if mask[ij]!=value # found the boundary
                out[IJ] = true
                break
            end
        end
    end
    return out
end

"""
    signed_distance(p::Point2{T}, poly::AbstractVector{Point2{T}}) signed_distance(p::Point2, poly::AbstractVector{<:Point2})where {T}

Returns the distance of a point `p` to the polygon `poly`.
"""
function signed_distance(p::Point2, poly::AbstractVector{<:Point2})
    if poly[1]==poly[end]
        poly = poly[1:end-1]
    end
    d = dot(p - poly[1], p - poly[1])
    s = 1.0
    j = length(poly)
    for i in eachindex(poly)
        e = poly[j] - poly[i]
        w = p - poly[i]
        b = w - e .* clamp(dot(w, e) / dot(e, e), 0.0, 1.0)
        d = min(d, dot(b, b))
        c = p[2] >= poly[i][2], p[2] < poly[j][2], e[1] * w[2] > e[2] * w[1]
        if all(c) || all(.!c)
            s = -s
        end
        j = i
    end
    return s * sqrt(d)
end


"""
    catchment_sinks(ctch_polygons::Vector{Vector{T}} where T<:Point2, routing_mask, x, y, sinks=get_boudary_cells(routing_mask, true);
                    check=true, threshold_dist=Inf,
                    D9_iters=4, D25_iter=2)
    catchment_sinks(ctch_polygons::AbstractVector, routing_mask, x, y, sinks=get_boudary_cells(routing_mask, true);
                    kws...)

Calculate the catchment sink cells for basins given by polygons.  The
output is intended to be fed to waterflows_subglacial via its `ctch_sinks` argument.
This is inteded to mark the cells of grounding line of a catchment as sinks.

Correctness check can be done with the `check_catchment_sinks` function, which is
done by default.

Input:
- ctch_polygons::Vector{AbstractVector{<:Point2}} -- a Vector of catchment polygons, presumably
                                             sourced from the internet.
                                             Need to be non-overlapping
- routing_mask - cells where water is routed set to true; assumed that boundary are all sinks
- x,y -- coordinate vectors
- sinks -- sink cells, by default calculated as all cells on the border of the routing_mask

KW-args
- check=true  -- do some checks (takes a bit of time)
- threshold_dist=Inf -- probably set to something lower

Output
- ctch_sinks -- a vector


Algo:
- get all boundary cells
  (make three vectors: cartesian-index, closest catchment, distance)
- calc signed_distance for all pixels for all catchments
- assign groundingline cells to catchment to which it is closest
- do some cleanup...
"""
catchment_sinks(ctch_polygons::AbstractVector, routing_mask, x, y, sinks=findall(get_boudary_cells(routing_mask, true)); kws...) =
    catchment_sinks([Point2.(cp) for cp in ctch_polygons], routing_mask, x, y, sinks; kws...)
function catchment_sinks(ctch_polygons::Vector{Vector{T}} where T<:Point2, routing_mask, x, y,
                         sinks=findall(get_boudary_cells(routing_mask, true));
                         check=true, threshold_dist=Inf,
                         D9_iters=4, D25_iter=2)
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

    ## Do cleanup
    # fix stray single gl-points: first do mode-filter in 9-neighborhood then 25-hood
    for _=1:D9_iters
        for s in sinks
            pts = [ctch_sinks_map[s]] # that way it gets picked up by `mode` in case of ties
            for ij in WWF.iterate_D9(s, ctch_sinks_map)
                ij == s && continue
                if ctch_sinks_map[ij]!=0
                    push!(pts, ctch_sinks_map[ij])
                end
            end
            ctch_sinks_map[s] = mode(pts)
        end
    end
    for _=1:D25_iter # too many iterations make things worse, 2 is optimal for IMBIE catchments in Antarctica
        for s in sinks
            pts = [ctch_sinks_map[s]] # that way it gets picked up by `mode` in case of ties
            for ij in iterate_D25(s, ctch_sinks_map)
                ij == s && continue
                if ctch_sinks_map[ij]!=0
                    push!(pts, ctch_sinks_map[ij])
                end
            end
            ctch_sinks_map[s] = mode(pts)
        end
    end

    ## Checks that
    check && check_catchment_sinks(ctch_sinks_map, ctch_polygons, sinks)

    ## turn into ctch_sinks which can be fed to waterflows_subglacial
    ctch_sinks = [findall(ctch_sinks_map.==i) for i = 1:length(ctch_polygons)]

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
