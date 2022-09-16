module WhereTheWaterFlows

using StaticArrays, Requires

export waterflows

const I11 = CartesianIndex(1,1)
const I22 = CartesianIndex(2,2)
const NOFLOW = 5     # direction number indicating no flow.  Use a constant to better keep track.
const NOFLOWer = 10  # direction number indicating no flow as well as no flow into this cell.

"""
Direction numbers.  E.g. dirnums[1,1] will return the number
corresponding to the direction top-left.

Note, I use the conversion that the x-axis corresponds to the row of
the matrix and the y-axis the columns.  To print them in this "normal"
coordinate system use `showme`
"""
const dirnums = SMatrix{3,3}(reverse([ 7      8 9
                                       4 NOFLOW 6
                                       1      2 3]',
                                     dims=2))

"Translation from dirnums to CartesianIndex"
const cartesian = SMatrix{3,3}(reverse(permutedims([CartesianIndex(-1,1)  CartesianIndex(0,1)  CartesianIndex(1,1)
                                                    CartesianIndex(-1,0)  CartesianIndex(0,0)  CartesianIndex(1,0)
                                                    CartesianIndex(-1,-1) CartesianIndex(0,-1) CartesianIndex(1,-1)]),
                                       dims=2))

"Adjustment factor for gradient in diagonal vs non-diagonal"
diagonal_fac(J::CartesianIndex) = (J.I[1]==0 || J.I[2]==0) ?  1.0 : 1/sqrt(2)


"""
Show an array with its axes oriented such that they correspond
to x and y direction.
"""
showme(ar) = (display(reverse(ar',dims=1)); println(" "))

"Translate a CartesianIndex, giving the offset, into a direction number"
ind2dir(ind::CartesianIndex) = dirnums[ind + I22]

"Translate a D8 direction number into a CartesianIndex (i.e. a flow vector). Maps dir==10 to no-flow also."
dir2ind(dir) = dir==10 ? CartesianIndex(0,0) : cartesian[dir]

"Translate a D8 direction number into a 2D vector."
dir2vec(dir) = [dir2ind(dir).I...]

"""
Tests whether a cell `J` with flowdir `dirJ` flows into cell `I`.
"""
flowsinto(J::CartesianIndex, dirJ::Integer, I::CartesianIndex) = ind2dir(I-J) == dirJ

"Return CartesianIndices corresponding to the 8 neighbors and the point itself"
iterate_D9(I, Iend, I1=I11) = max(I1, I-I1):min(Iend, I+I1)
function iterate_D9(I, ar::AbstractMatrix)
    R = CartesianIndices(size(ar))
    I1, Iend = first(R), last(R)
    return max(I1, I-I1):min(Iend, I+I1)
end

"""
    on_outer_boundary(ar, I::CartesianIndex)

Check whether on outer boundary of ar
"""
function on_outer_boundary(ar, I::CartesianIndex)
    i,j = I.I
    iend, jend = size(ar)
    if i==1 || i==iend || j==1 || j==jend
        return true
    else
        return false
    end
end

"""
    d8dir_feature(dem, bnd_as_pits)

D8 directions of a DEM and drainage features.

Elevations with NaN map to dir==NOFLOW, i.e. just like pits.
However, they are not treated in the pits-filling function `drainpits`.

The argument `bnd_as_pits` determines whether neighboring cells
flow out of the domain or into NaN-cells (`bnd_as_pits==true`)
or not.

Return
- dir  - direction, encoded as `dirnums`
- nout - number of outflow cells of a cell (0 or 1)
- nin  - number of inflow cells of a cell (0-8)
- pits - location of pits as a `Vector{CartesianIndex{2}}` (sorted)
- flowdir_extra_output -- nothing (not used by this function)
"""
function d8dir_feature(dem, bnd_as_pits)
    # outputs
    dir = fill!(similar(dem, Int8), 0)
    #nout = falses(size(dem))
    nout = fill!(similar(dem, Bool), false)
    nin = fill!(similar(dem, Int8), 0)
    pits = CartesianIndex{2}[]

    R = CartesianIndices(size(dem))
    Iend = last(R)

    # get dir for all points
    for I in R
        # make pits on boundary if bnd_as_pits is set
        if bnd_as_pits && on_outer_boundary(dem,I)
            # make it a pit
            dir[I] = NOFLOW
            continue
        end

        ele = dem[I]
        delta_ele = 0.0 # keeps track of biggest elevation change
        dir_ = NOFLOW
        if isnan(ele)
            # just mark as NOFLOW
        else
            for J in iterate_D9(I, Iend)
                I==J && continue
                ele2 = dem[J]
                if isnan(ele2)
                    if bnd_as_pits
                        # flow into first found NaN-cell
                        dir_ = ind2dir(J-I)
                        break
                    else
                        # ignore NaN-Cell
                        continue
                    end
                end
                delta_ele2 = (ele2 - ele) * diagonal_fac(J-I)
                if delta_ele2 < delta_ele
                    # lower elevation found, adjust dir
                    delta_ele = delta_ele2
                    dir_ = ind2dir(J-I)
                end
            end
        end
        dir[I] = dir_
    end
    # flow features
    # (not thread-safe because of push!)
    for I in R
        for J in iterate_D9(I, Iend)
            J==I && continue
            nin[I] += flowsinto(J, dir[J], I)
        end
        if dir[I]==NOFLOW
            nout[I] = false
            if !isnan(dem[I])
                push!(pits, I)
            elseif nin[I]>0
                # mark NaN points only as pits if something is flowing into them
                push!(pits, I)
            end
        else
            nout[I] = true
        end
    end

    return dir, nout, nin, pits, dem, nothing
end

"""
    waterflows(dem, cellarea=cellarea=fill!(similar(dem),1), flowdir_fn=d8dir_feature;
               calc_streamlength=true, drain_pits=true, bnd_as_pits=false)

Does the water flow routing according the D8 algorithm.  Locations of the `dem`
with `NaN`-value are ignored.

args:
- dem -- the DEM (or hydro-potential); array
- cellarea=cellarea=fill!(similar(dem),1) -- the source per cell, defaults to 1.  If using physical units
                              then use a volumetric flux per cell, e.g. m3/s.
- flowdir_fn=d8dir_feature -- the routing function.  Defaults to the built-in `d8dir_feature`
                              function but could be customized

kwargs:
- drain_pits -- whether to route through pits (true)
- bnd_as_pits (true) -- whether the domain boundary and NaNs should be pits,
                 i.e. adjacent cells can drain into them,
                 or whether to ignore them.

TODO: if bnd_as_pits routes water along boundary edges first. This would probably
substantially reduce the number of catchments, as currently every boundary point
is a pit and thus a catchment (if bnd_as_pits==true).

Returns
- area -- upslope area
- stream-length -- length of stream to the farthest source
- dir -- flow direction at each location
- nout -- whether the point has outlflow.  I.e. nout[I]==0 --> I is a pit
- nin -- number of inflow cells
- pits -- location of pits as Vector{CartesianIndex{2}}
- c -- catchment map
- bnds -- boundaries between catchments.  The boundary to the exterior/NaNs is not in here.
- flowdir_extra_output -- extra output of the flowdir_fn, which is `nothing` for the default
"""
function waterflows(dem, cellarea=fill!(similar(dem),1), flowdir_fn=d8dir_feature;
                    drain_pits=true, bnd_as_pits=true)
    dir, nout, nin, pits, dem4drainpits, flowdir_extra_output = flowdir_fn(dem, bnd_as_pits)
    area, slen, c = flowrouting_catchments(dir, pits, cellarea)
    bnds = make_boundaries(c, 1:length(pits))
    if drain_pits
        drainpits!(dir, nin, nout, pits, c, bnds, dem4drainpits)
        area, slen, c = flowrouting_catchments(dir, pits, cellarea)
    end
    #area[isnan.(dem)] .= NaN
    return area, slen, dir, nout, nin, pits, c, bnds, flowdir_extra_output
end

# TODO: implement this instead of waterflows
# """
#     flowrouting(dir, nin, cellarea; maxiter=max(size(dir)...)*2, calc_streamlength=true)

# Do the actual routing

# Input:
# - `dir` direction field
# - `nin` number of inputs into cell
# - `cellarea` water input per cell

# KW:
# - `maxiter=max(size(dir)...)*2` -- may need to be increased for very sinous drainage paths
# - `calc_streamlength=true` -- calculate stream length.  Will slow-down calculation somewhat

# Return
# - upstream area
# - stream-length (or something irrelevant if calc_streamlength==false)
# """

"""
    flowrouting_catchments(dir, pits, cellarea)

Recursively calculate flow-routing and catchments from
- dir - direction field
- pits - pit coordinates
- cellarea - water input

Return:
- upstream area
- stream length
- catchments Matrix{Int}.  Value==0 corresponds to NaNs in the DEM
  which are not pits (i.e. where no water flows into).

Note: this function may cause a stackoverflow on very big catchments.
"""
function flowrouting_catchments(dir, pits, cellarea)
    c = fill!(similar(dir, Int), 0)
    area = fill!(similar(dir,Float64), 0)
    slen = fill!(similar(dir, Int), 0)
    np = length(pits)
    # recursively traverse the drainage tree in up-flow direction,
    # starting at all pits
    Threads.@threads for color = 1:length(pits)
        pit = pits[color]
        _flowrouting_catchments!(area, slen, c, dir, cellarea, color, pit)
    end
    return area, slen, c
end
# modifies c and area
function _flowrouting_catchments!(area, len, c, dir, cellarea, color, ij)
    c[ij] = color

    # proc upstream points
    uparea = 0.0
    slen = 0
    for IJ in iterate_D9(ij, c)
        if ij==IJ
            uparea += cellarea[ij]
            slen = max(slen, 1)
        elseif flowsinto(IJ, dir[IJ], ij)
            uparea_, slen_ = _flowrouting_catchments!(area, len, c, dir, cellarea, color, IJ)
            uparea += uparea_
            slen = max(slen, slen_+1)
        end
    end
    area[ij] = uparea
    len[ij] = slen
    return uparea, slen
end

"""
    make_boundaries(catchments, colors)

Make vectors of boundary points.  Note that points along the
edge of the domain as well as points bordering NaN-cells with no
inflow do not count as boundaries.
"""
function make_boundaries(catchments, colors)
    bnds = [CartesianIndex{2}[] for c in colors]
    for R in CartesianIndices(size(catchments))
        c = catchments[R]
        c==0 && continue # don't find boundaries for c==0 (NaNs with no inflow)
        bnd = bnds[c]
        for I in iterate_D9(R, catchments)
            co = catchments[I]
            if co!=c && co!=0
                push!(bnd, R)
                break
            end
        end
    end
    bnds
end

"""
Checks that all points in `bnds[color]` are indeed on the
boundary; removes them otherwise.

TODO: this is one of the performance bottlenecks.
"""
function _prune_boundary!(del, bnds, catchments::AbstractMatrix, color, colormap)
    empty!(del)
    # loop over all boundary cells
    @inbounds for (i,P) in enumerate(bnds[color])
        keep = false
        # check cells around it
        for PP in iterate_D9(P, catchments)
            # if one is of different color, keep it.
            co = catchments[PP]
            if co!=0
                co = colormap[co]
                if co!=color
                    keep = true
                    break
                end
            end
        end
        if !keep
            push!(del, i)
        end
    end
    deleteat!(bnds[color], del)
    return nothing
end

"""
     drainpits!(dir, nin, nout, pits, c, bnds, dem;
               maxiter=100)

Return an updated direction field which drains (interior) pits.
This is done by reversing the flow connecting the lowest point on the catchment boundary
to the pit for each catchment.

There needs to be a decision on how to treat the outer boundary and boundaries
to NaN-cells.  Possibilities:
- do not treat boundary cells as pits but do treat pits at the boundary as terminal.
- treat boundary cells as pits, i.e. flow reaching such a cell will vanish.  Again
  such cells would be terminal.  This is currently done

TODO: What to do if there are no terminal pits in the DEM?  Fill to the uppermost pit?  Or take
the lowermost as terminal? -> currently it picks a random one

Returns new dir, nin, nout, pits (sorted), c, bnds

TODO: this is the performance bottleneck.
"""
function drainpits!(dir, nin, nout, pits, c, bnds, dem)
    maxiter = 100
    Iend = CartesianIndex(size(dem))

    pits_to_keep = trues(length(pits))

    # Table which translates the old color to the new color.
    # Note that only the currently processed color might change.
    # Initialize to the id-map (0 is used for off-points)
    colormap = collect(1:length(pits))

    colormap_inv = [[i] for i=1:length(pits)] # note that the total length

    no_drainage_across_boundary = false

    # work array as otherwise much time is spent allocating
    # (if below loop is ever multi-threaded, this will need one per thread)
    del_workarray  = Int[]
    sizehint!(del_workarray, maximum(length.(bnds)))

    # iterate until all interior pits are removed (not really needed, I think)
    @inbounds for i=1:maxiter
        n_removed = 0
        for (color, P) in enumerate(pits)
            # Already removed pit, skip
            P==CartesianIndex(-1,-1) && continue
            # Don't process pits on the DEM boundary, because there water disappears.
            (P.I[1]==1 || P.I[2]==1 || P.I[1]==size(dir,1) || P.I[2]==size(dir,2)) && continue
            # Don't process pits which are NaNs, because there water disappears.
            # See also bnd_as_pits.
            isnan(dem[P]) && continue
            # If there are no more boundaries left, stop.  This should only occur when
            # bnd_as_pits==false and when there are only interior pits.
            if all((isempty(b) for b in bnds))
                no_drainage_across_boundary = true # set flag to warn later
                break
            end
            # If there are no boundaries for this point, go to next point
            isempty(bnds[color]) && continue

            # Find point on catchment boundary with minimum elevation
            # Avoid, if possible, picking a NOFLOWer point.
            minn = typemax(eltype(dem))
            Imin = CartesianIndex(-1,-1)
            for I in bnds[color]
                if dem[I] < minn && dir[I]!=NOFLOWer
                    minn = dem[I]
                    Imin = I
                end
            end
            if Imin == CartesianIndex(-1,-1) # no NOFLOWer cells found
                for I in bnds[color]
                    if dem[I] < minn
                        minn = dem[I]
                        Imin = I
                    end
                end
            end

            @assert colormap[c[Imin]]==c[P]==color # Something is amiss if the found minimum is the pit!
            # make the outflow and find the next catchment
            min_ = Inf
            target = Imin
            for J in iterate_D9(Imin, dem)
                Imin==J && continue # do not look at the point itself
                cc = c[J]
                cc = cc==0 ? 0 : colormap[cc]
                cc == color && continue # J is in Imin's catchment
                if dem[J] < min_
                    min_ = dem[J]
                    target = J
                end
            end
            if target==Imin
                # this means that point is on a NaN boundary -> don't process
                continue
            end

            # Reverse directions on path going from Imin to P
            P1 = Imin
            P2 = dir2ind(dir[P1]) + P1
            # first, do the flow across the boundary
            _flow_from_to!(Imin, target, dir, nin, nout)
            while P1!=P
                # get down-downstream point (do ahead lookup as a undisturbed dir is needed)
                P3 = dir2ind(dir[P2]) + P2
                # reverse
                _flow_from_to!(P2, P1, dir, nin, nout)

                P1 = P2
                P2 = P3
            end

            # remove from list of pits
            pits_to_keep[color] = false
            pits[color] = CartesianIndex(-1,-1)

            # update colormap and boundaries
            othercolor = colormap[c[target]]
            ## this is slow
            # Threads.@threads for i in eachindex(colormap)
            #     if color==colormap[i]
            #         colormap[i] = othercolor
            #     end
            # end
            ## thus do this instead:
            append!(colormap_inv[othercolor], colormap_inv[color])
            empty!(colormap_inv[color])
            for color in colormap_inv[othercolor]
                colormap[color] = othercolor
            end

            append!(bnds[othercolor], bnds[color])
            empty!(bnds[color])
            _prune_boundary!(del_workarray, bnds, c, othercolor, colormap)
            n_removed +=1
        end
        n_removed==0 && break
        if i==maxiter
            error("Maximum number of iterations reached in `drainpits`: $i")
        end
    end

    if no_drainage_across_boundary
        @warn """Water cannot leave this DEM!  Instead it drains into a random interior pit.
                 Consider setting `bnd_as_pits=true`."""
    end

    # make new colors consecutive
    newcolors = unique(colormap)
    d = Dict((c=>i for (i,c) in enumerate(newcolors)))
    colormap2 = [d[c] for c in colormap]

    # remove removed pits and bnds
    deleteat!(pits, .!pits_to_keep)
    deleteat!(bnds, .!pits_to_keep)

    # # update catchments -> done with `flowrouting_catchments`
    return nothing
end

"""
Update dir, nin, and nout such that flow at P1 is now
from P1 to P2.

It can potentially modify `dir, nin, nout` at three locations
- P1: dir, nout, nin
- P2: nin
  - if allow_reversion==true, then P2's dir, nout can also be modified.
- P3 (previous receiver cell of P1): nin

Note that if P1 and P2 lie in the same catchment, then P3 is also in that catchment.
"""
function _flow_from_to!(P1, P2, dir, nin, nout, allow_reversion=false)
    # already right
    ind2dir(P2-P1)==dir[P1] && return nothing

    # otherwise update
    P3 = P1 + dir2ind(dir[P1]) # previous receiver cell
    if P3!=P1
        nin[P3] -= 1
    end
    dir[P1] = ind2dir(P2-P1)
    nin[P2] += 1
    nout[P1] = true # definitely not a pit anymore

    # if flow from P2 was into P1, then make P2 a pit
    # (otherwise dir will be inconsistent)
    if (dir2ind(dir[P2]) + P2 == P1)
        if allow_reversion
            dir[P2] = ind2dir(CartesianIndex(0,0))
            nout[P2] = false
            nin[P1] -= 1
        else
            error("Flow direction in P2 would need to be reversed but not allowed (set option `allow_reversion`).")
        end
    end
    return nothing
end


## Plotting
function __init__()
    @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" include("plotting.jl")
end

## Post processing
include("postproc.jl")

end # module
