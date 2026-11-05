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
    Base.require_one_based_indexing(A, window, weights)
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
    get_boundary_cells(mask, value)

Get all the cells of mask of value `value` which border
on cells of different value.
"""
function get_boundary_cells(mask, value)
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

