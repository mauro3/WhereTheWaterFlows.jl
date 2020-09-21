module WhereTheWaterFlows

using StaticArrays, Requires

export waterflows, fill_dem, catchments,
    plotarea

const I11 = CartesianIndex(1,1)
const I22 = CartesianIndex(2,2)
#const NaNflow = -99
const NOFLOW = 5 # direction number indicating no flow.  Use a constant to better keep track.

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
# "As dirnums but with no-flow encoded as NaNflow."
# const dirnums_nan = SMatrix{3,3}(reverse([ 7       8 9
#                                            4 NaNflow 6
#                                            1       2 3]',
#                                          dims=2))

"Translation from dirnums to CartesianIndex"
const cartesian = SMatrix{3,3}(reverse(permutedims([CartesianIndex(-1,1)  CartesianIndex(0,1)  CartesianIndex(1,1)
                                                    CartesianIndex(-1,0)  CartesianIndex(0,0)  CartesianIndex(1,0)
                                                    CartesianIndex(-1,-1) CartesianIndex(0,-1) CartesianIndex(1,-1)]),
                                       dims=2))
"""
Show an array with its axes oriented such that they correspond
to x and y direction.
"""
showme(ar) = (display(reverse(ar',dims=1)); println(" "))

"Translate a CartesianIndex, giving the offset, into a direction number"
ind2dir(ind::CartesianIndex) = dirnums[ind + I22]
# ind2dir_nan(ind::CartesianIndex) = dirnums_nan[ind + I22]

"Translate a D8 direction number into a CartesianIndex"
dir2ind(dir) = cartesian[dir]
# "Translate a D8 direction number (with NaNflow at center) into a CartesianIndex"
# dir2ind_nan(dir) = dir==0 ? CartesianIndex(0,0) : dir2ind(dir)

"Translate a D8 direction number into a 2D vector."
dir2vec(dir) = [dir2ind(dir).I...]

"""
Tests whether a cell `J` with flowdir `dirJ` flows into cell `I`.
"""
flowsinto(J::CartesianIndex, dirJ, I::CartesianIndex) = ind2dir(I-J) == dirJ

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
"""
function d8dir_feature(dem, bnd_as_pits)
    # outputs
    diro = zeros(Int8, size(dem))
    nout = falses(size(dem))
    nin = zeros(Int8, size(dem))
    pits = CartesianIndex{2}[]

    R = CartesianIndices(size(dem))
    Iend = last(R)

    # get dir for all points
    for I in R
        # make pits on boundary if bnd_as_pits is set
        if bnd_as_pits && on_outer_boundary(dem,I)
            # make it a pit
            diro[I] = NOFLOW
            continue
        end

        ele = dem[I] # keeps track of lowest elevation of all 9 cells
        dir = NOFLOW
        if isnan(ele)
            # just mark as NOFLOW
        else
            for J in iterate_D9(I, Iend)
                I==J && continue
                ele2 = dem[J]
                if isnan(ele2)
                    if bnd_as_pits
                        # flow into first found NaN-cell
                        dir = ind2dir(J-I)
                        break
                    else
                        # ignore NaN-Cell
                        continue
                    end
                elseif ele > ele2
                    # lower elevation found, adjust dir
                    ele = ele2
                    dir = ind2dir(J-I)
                end
            end
        end
        diro[I] = dir
    end
    # flow features
    for I in R
        for J in iterate_D9(I, Iend)
            J==I && continue
            nin[I] += flowsinto(J, diro[J], I)
        end
        if diro[I]==NOFLOW
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

    return diro, nout, nin, pits
end

"""
    waterflows(dem, cellarea=ones(size(dem));
               maxiter=1000, calc_streamlength=true, drain_pits=true, bnd_as_pits=false)

Does the water flow routing according the D8 algorithm.

kwargs:
- drain_pits -- whether to route through pits (true)
- maxiter -- maximum iterations of the algorithm (1000)
- calc_streamlength -- whether to calculate stream length (true)
- bnd_as_pits (true) -- whether the domain boundary and NaNs should be pits,
                 i.e. adjacent cells can drain into them,
                 or whether to ignore them.

TODO: if bnd_as_pits rout water along boundary edges first. This would probably
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
"""
function waterflows(dem, cellarea=ones(size(dem));
                    maxiter=1000, calc_streamlength=true, drain_pits=true, bnd_as_pits=true)
    area, slen, dir, nout, nin, pits = _waterflows(d8dir_feature(dem, bnd_as_pits)..., cellarea; maxiter=maxiter, calc_streamlength=calc_streamlength)
    c, bnds = catchments(dir, pits, bnd_as_pits)
    if drain_pits
        dir, nin, nout, pits, c, bnds = drainpits(dem, dir, nin, nout, pits, (c,bnds))
        area, slen = _waterflows(dir, nout, nin, pits, cellarea; maxiter=maxiter, calc_streamlength=calc_streamlength)
    end
    #area[isnan.(dem)] .= NaN
    return area, slen, dir, nout, nin, pits, c, bnds
end
# this function does the actual routing
function _waterflows(dir, nout, nin, pits, cellarea=ones(size(dir)); maxiter=1000, calc_streamlength=true)
    area = copy(cellarea)
    tmp = convert(Matrix{Int}, nin)
    tmp2 = calc_streamlength ? copy(tmp) : tmp


    for counter = 1:maxiter
        n = 0
        for R in CartesianIndices(size(dir))
            # Consider only points with less than `counter` upstream points
            if tmp[R]==0
                tmp[R] = -counter # done with it
                if calc_streamlength
                    tmp2[R] = -counter # done with it
                end
                d = dir[R]
                if d!=NOFLOW
                    receiver = R + dir2ind(d)
                    area[receiver] += area[R]
                    tmp2[receiver] -= 1
                    tmp2[receiver] < 0 && error("This should not happen!")
                end
                n +=1
            end
        end
        if n==0
            #@show counter
            break
        end
        if counter==maxiter
            error("Maximum number of iterations reached in `_waterflows`: $counter")
        end
        copyto!(tmp, tmp2)
    end
    return area, -tmp, dir, nout, nin, pits
end

"""
    catchments(dir, pits, bnd_as_pits)

Calculate catchments from
- dir
- pits

Return: catchments Matrix{Int}.  Value==0 corresponds to NaNs in the DEM
which are not pits (i.e. where no water flows into).
"""
function catchments(dir, pits, bnd_as_pits)
    c = zeros(Int, size(dir))
    np = length(pits)
    # recursively traverse the drainage tree in up-flow direction,
    # starting at all pits
    for (n, pit) in enumerate(pits)
        _catchments!(n, c, dir, pit)
    end

    return c, make_boundaries(c, 1:np)
end
function _catchments!(n, c, dir, ij)
    c[ij] = n
    # proc upstream points
    for IJ in iterate_D9(ij, c)
        ij==IJ && continue
        if flowsinto(IJ, dir[IJ], ij)
            _catchments!(n, c, dir, IJ)
        end
    end
end

"""
    make_boundaries(catchments, colors, bnd_as_pits)

Make vectors of boundary points.  Note that points along the
edge of the domain as well as points bordering NaN-cells with no
inflow do not count as boundaries.

TODO: this is brute force...
"""
function make_boundaries(catchments, colors)
    bnds = [CartesianIndex[] for c in colors]
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
    end#
    bnds
end

"""
Checks that all points in `bnds[color]` are indeed on the
boundary; removes them otherwise.
"""
function _prune_boundary!(bnds, catchments::Matrix, color)
    del = Int[]
    # loop over all boundary cells
    for (i,P) in enumerate(bnds[color])
        keep = false
        # check cells around it
        for PP in iterate_D9(P, catchments)
            # if one is of different color, keep it.
            co = catchments[PP]
            if co!=color && co!=0
                keep = true
                break
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
     drainpits(dem, dir, nin, nout, pits, (c, bnds)=catchments(dir, pits);
               maxiter=100)

Return an update direction field which drains (interior) pits.
This is done by reversing the flow connecting the lowest boundary point
to the pits.

There needs to be a decision on how to treat the outer boundary and boundaries
to NaN-cells.  Possibilities:
- do not treat boundary cells as pits but do treat pits at the boundary as terminal.
- treat boundary cells as pits, i.e. flow reaching such a cell will vanish.  Again
  such cells would be terminal.

What to do if there are no terminal pits in the DEM?  Fill to the uppermost pit?  Or take
the lowermost as terminal?

Returns new dir, nin, nout, pits (sorted), c.

TODO: this is the performance bottleneck.
"""
function drainpits(dem, dir, nin, nout, pits, (c, bnds)=catchments(dir, pits);
                   maxiter=100)
    dir_ = copy(dir)
    dir2 = copy(dir)
    nin_ = copy(nin)
    nout_ = copy(nout)
    pits_ = copy(pits)
    c_ = copy(c)
    bnds_ = deepcopy(bnds)

    Iend = CartesianIndex(size(dem))

    pits_to_keep = trues(length(pits_))

    no_drainage_across_boundary = false
    # iterate until all interior pits are removed
    for i=1:maxiter
        n_removed = 0
        for (color, P) in enumerate(pits_)

            # Already removed pit, skip
            P==CartesianIndex(-1,-1) && continue
            # Don't process pits on the DEM boundary, because there water disappears.
            (P.I[1]==1 || P.I[2]==1 || P.I[1]==size(dir,1) || P.I[2]==size(dir,2)) && continue
            # Don't process pits which are NaNs, because there water disappears.
            # See also bnd_as_pits.
            isnan(dem[P]) && continue
            # If there are no more boundaries left, stop.  This should only occur when
            # bnd_as_pits==false and when there are only interior pits.
            if all(isempty.(bnds_))
                no_drainage_across_boundary = true # set flag to warn later
                break
            end
            # If there are no boundaries for this point, go to next point
            isempty(bnds_[color]) && continue

            ## Debug plotting
            # cls = [c=>sum(c_.==c) for c=1:6]
            # @show color, size(bnds_[color])
            # @show cls
            # imshow(Array(c_'), origin="lower")
            # if color==1
            #     colorbar()
            # end
            # readline(stdin)
            # if color==6
            #     imshow(Array(c_'), origin="lower")
            # end

            # find point on catchment boundary with minimum elevation
            Imin = bnds_[color][findmin(getindex.(Ref(dem), bnds_[color]))[2]]
            @assert c_[Imin]==c_[P] # Something is amiss if the found minimum is the pit!
            # make the outflow and find the next catchment
            min_ = Inf
            target = Imin
            for J in iterate_D9(Imin, dem)
                Imin==J && continue # do not look at the point itself
                c_[J] == color && continue # J is in Imin's catchment
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
            P2 = dir2ind(dir_[P1]) + P1
            # first, do the flow across the boundary
            _flow_from_to!(Imin, target, dir_, nin_, nout_)
            while P1!=P
                # get down-downstream point (do ahead lookup as a undisturbed dir_ is needed)
                P3 = dir2ind(dir_[P2]) + P2
                # reverse
                _flow_from_to!(P2, P1, dir_, nin_, nout_)

                P1 = P2
                P2 = P3
            end

            # remove from list of pits
            pits_to_keep[color] = false
            pits_[color] = CartesianIndex(-1,-1)

            # update catchments
            othercolor = c_[target]
            c_[c_.==color] .= othercolor
            append!(bnds_[othercolor], bnds_[color])
            empty!(bnds_[color])
            _prune_boundary!(bnds_, c_, othercolor)
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

    # remove removed pits
    pits_ = pits_[pits_.!=Ref(CartesianIndex(-1, -1))]
    # redo colors
    cols = Dict([c=>i for (i,c) in enumerate(unique(c_))])
    for i in eachindex(c_)
        c_[i] = cols[c_[i]]
    end

    return dir_, nin_, nout_, sort(pits_), c_, bnds_
end

"""
Update dir, nin, and nout such that flow at P1 is now
from P1 to P2.

It can potentially modify `dir, nin, nout` at three locations
- P1: dir, nout, nin
- P2: nin
  - if allow_reversion==true, then P2's dir, nout can also be modified.
- P3 (previous receiver cell of P1): nin
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

"""
    fill_dem(dem, pits, dir)

Fill the pits (aka sinks) of a DEM (apply this after applying "drainpits",
which is the default). Returns the filled DEM.

Note, this is not needed as pre-processing step to use `upstream` area.

This uses a tree traversal to fill the DEM. It does it depth-first (as it
is easier) which may lead to a stack overflow on a large DEM.
"""
function fill_dem(dem, pits, dir; small=0.0)
    dem = copy(dem)
    for pit in pits
        npts = _fill_ij!(0.0, dem, pit, dir, small)
    end
    return dem
end

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

## Plotting
function __init__()
    @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" include("plotting.jl")
end

end # module
