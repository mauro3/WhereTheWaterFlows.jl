module WhereTheWaterFlows

export waterflows, fill_dem, catchments,
    plotarea

# Note, I use the conversion (as in VAWTools) that the x-axis corresponds
# to the row of the matrix and the y-axis the columns.  To print them in this "normal"
# coordinate system use `showme`
using PyPlot, StaticArrays

"""
Direction numbers.  E.g. dirnums[1,1] will return the number
corresponding to the direction top left.

Note, I use the conversion that the x-axis corresponds to the row of
the matrix and the y-axis the columns.  To print them in this "normal"
coordinate system use `showme`
"""
const dirnums = SMatrix{3,3}(reverse([ 7 8 9
                                       4 5 6
                                       1 2 3]', dims=2))
"As dirnums but with dirnums..."
const dirnums_0 = SMatrix{3,3}(reverse([ 7 8 9
                                         4 0 6
                                         1 2 3]', dims=2))
"""
Show an array with its axes oriented such that they correspond
to x and y direction.
"""
showme(ar) = (display(reverse(ar',dims=1)); println(" "))

const I11 = CartesianIndex(1,1)
const I22 = CartesianIndex(2,2)

"""
Tests whether a cell `J` with flowdir `dirJ` flows into cell `I`.
"""
flowsinto(J::CartesianIndex, dirJ, I::CartesianIndex) = dirnums_0[I-J+I22] == dirJ

"Translate a D8 direction number into a CartesianIndex"
dir2ind(dir) = findfirst(x->x==dir, dirnums) - I22
"Translate a D8 direction number into a 2D vector."
dir2vec(dir) = [dir2ind(dir).I...]
"Translate a CartesianIndex, giving the offset, into a direction number"
ind2dir(ind::CartesianIndex) = dirnums[ind + I22]

"Return CartesianIndices corresponding to the 8 neighbors and the point itself"
iterate_D9(I, Iend, I1=I11) = max(I1, I-I1):min(Iend, I+I1)
function iterate_D9(I, ar::AbstractMatrix)
    R = CartesianIndices(size(ar))
    I1, Iend = first(R), last(R)
    return max(I1, I-I1):min(Iend, I+I1)
end

"""
    d8dir_feature(dem)

D8 directions of a DEM and drainage features.

Return
- dir - direction, see `dirnums`
- nout - number of  outflow cells of a cell (0 or 1)
- nin - number of inflow cells of a cell (0-8)
"""
function d8dir_feature(dem)
    # TODO this does not do the "right" thing at the edges.  Although, maybe that is just fine.

    diro = zeros(Int8, size(dem))
    nout = falses(size(dem))
    nin = zeros(Int8, size(dem))
    sx,sy = size(dem)
    # https://julialang.org/blog/2016/02/iteration
    # https://discourse.julialang.org/t/psa-replacement-of-ind2sub-sub2ind-in-julia-0-7/14666
    R = CartesianIndices(size(dem))
    Iend = last(R)
    # inner points
    for I in CartesianIndices((1:sx, 1:size(dem,2)))
        ele = dem[I]
        dir = 5
        for J in iterate_D9(I, Iend)
            ele2 = dem[J]
            if ele > ele2
                ele = ele2
                dir = ind2dir(J-I) # equivalent to dir = i
                # @assert ind2dir(J-I)==i
            end
        end
        diro[I] = dir
    end
    # flow features
    for I in R
        nout[I] = diro[I]==5 ? false : true
        for J in iterate_D9(I, Iend)
            nin[I] += flowsinto(J, diro[J], I)
        end
    end

    pits = findall(nout.==false)
    return diro, nout, nin, pits
end

## This does the right thing at the edges but is ugly; and probably buggy.
# function d8dir_feature1(dem)
#     diro = zeros(Int8, size(dem))
#     nout = falses(size(dem))
#     nin = zeros(Int8, size(dem))
#     sx,sy = size(dem)
#     # https://julialang.org/blog/2016/02/iteration
#     # https://discourse.julialang.org/t/psa-replacement-of-ind2sub-sub2ind-in-julia-0-7/14666/8
#     R = CartesianIndices(size(dem))
#     I1, Iend = first(R), last(R)
#     # inner points
#     for I in CartesianIndices((2:sx-1, 2:size(dem,2)-1))
#         ele = dem[I]
#         dir = 5
#         for (i,J) in enumerate(iterate_D9(I, Iend, I1))
#             ele2 = dem[J]
#             if ele > ele2
#                 ele = ele2
#                 dir = ind2dir(J-I) # equivalent to dir = i
#                 # @assert ind2dir(J-I)==i
#             end
#         end
#         diro[I] = dir
#     end
#     # along edges
#     for x=1:2
#         dirs = ((7,8,9), (1,2,3))[x]
#         iy = (1,sy)[x]
#         iyy = (2,sy-1)[x]
#         for ix=2:sx-1
#             ele = dem[ix,iy]
#             dir = 5
#             for (i, ixx) = zip(dirs, ix-1:ix+1)
#                 ele2 = dem[ixx,iyy]
#                 if ele > ele2
#                     ele = ele2
#                     dir = i
#                 end
#             end
#             diro[ix,iy] = dir
#         end
#     end

#     for y=1:2
#         dirs = ((3,6,9), (1,4,7))[y]
#         ix = (1,sx)[y]
#         ixx = (2,sx-1)[y]
#         for iy=2:sy-1
#             ele = dem[ix,iy]
#             dir = 5
#             for (i, iyy) = zip(dirs, iy-1:iy+1)
#                 ele2 = dem[ixx,iyy]
#                 if ele > ele2
#                     ele = ele2
#                     dir = i
#                 end
#             end
#             diro[ix,iy] = dir
#         end
#     end

#     # corners
#     diro[1,1] = dem[1,1]>dem[2,2] ? 9 : 5
#     diro[1,sy] = dem[1,sy]>dem[2,sy-1] ? 3 : 5
#     diro[sx,1] = dem[sx,1]>dem[sx-1,2] ? 7 : 5
#     diro[sx,sy] = dem[sx,sy]>dem[sx-1,sy-1] ? 1 : 5

#     # flow features
#     for I in R
#         nout[I] = diro[I]==5 ? false : true
#         for J in iterate_D9(I, Iend)
#             # if I==CartesianIndex(1,1)
#             #     @show J, diro[J], I
#             #     @show flowsinto(J, diro[J], I)
#             # end
#             nin[I] += flowsinto(J, diro[J], I)
#         end
#     end

#     return diro, nout, nin, findall(nout.==false)
# end


"""
    waterflows(dem, cellarea=ones(size(dem))
    waterflows(dir, nout, nin, ic_area=ones(size(dir)), calc_streamlength=true)

kwargs:
- fillpits -- whether to route through pits
- maxiter -- maximum iterations of the algorithm
- calc_streamlength -- whether to calculate stream length

Returns
- upslope area (for a D8 direction field)
- stream-length -- number of upstream cells for each cell
- dir -- flow direction at each location
- nout -- whether the point has outlflow.  I.e. nout[I]==0 --> I is a pit
- nin -- number of inflow cells
- pits -- location of pits as Vector{CartesianIndex{2}}
- c -- catchment map
- bnds -- boundaries between catchments
"""
function waterflows end
function waterflows(dem, cellarea=ones(size(dem));
                     maxitr=1000, calc_streamlength=true, fillpits=true)
    area, slen, dir, nout, nin, pits = _waterflows(d8dir_feature(dem)..., cellarea; maxitr=maxitr, calc_streamlength=calc_streamlength)
    c, bnds = catchments(dir, pits)
    if fillpits
        dir, nin, nout, pits, c, bnds = drainpits(dem, dir, nin, nout, pits, (c,bnds))
        area, slen = _waterflows(dir, nout, nin, pits)
    end
    return area, slen, dir, nout, nin, pits, c, bnds
end
function _waterflows(dir, nout, nin, pits, cellarea=ones(size(dir)); maxitr=1000, calc_streamlength=true)
    area = copy(cellarea)
    tmp = convert(Matrix{Int}, nin)
    tmp2 = calc_streamlength ? copy(tmp) : tmp
    for counter = 1:maxitr
        n = 0
        for R in CartesianIndices(size(dir))
            # Consider only points with less than `counter` upstream points
            if tmp[R]==0
                tmp[R] = -counter # done with it
                if calc_streamlength
                    tmp2[R] = -counter # done with it
                end
                d = dir[R]
                if d!=5
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
        copyto!(tmp, tmp2)
    end
    return area, -tmp, dir, nout, nin, pits
end

"""
    catchments(dir, pits)

Calculate catchments from
- dir
- pits
"""
function catchments(dir, pits)
    c = zeros(Int, size(dir))
    np = length(pits)
    # recursively traverse the drainage tree
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
Make vectors of boundary points.

TODO: this is brute force...
"""
function make_boundaries(catchments, colors)
    bnds = [CartesianIndex[] for c in colors]
    for R in CartesianIndices(size(catchments))
        c = catchments[R]
        bnd = bnds[c]
        for I in iterate_D9(R, catchments)
            if catchments[I]!=c
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
"""
function _prune_boundary!(bnds, catchments::Matrix, color)
    del = Int[]
    # loop over all boundary cells
    for (i,P) in enumerate(bnds[color])
        keep = false
        # check cells around it
        for PP in iterate_D9(P, catchments)
            # if one is of different color, keep it.
            if catchments[PP]!=color
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
Return an update direction field which drains pits.

Returns new dir, nin, nout, pits, c.

TODO: this is the performance bottleneck.
"""
function drainpits(dem, dir, nin, nout, pits, (c, bnds)=catchments(dir, pits);
                   maxitr=100)
    dir_ = copy(dir)
    dir2 = copy(dir_)
    nin_ = copy(nin)
    nout_ = copy(nout)
    pits_ = copy(pits)
    c_ = copy(c)
    bnds_ = deepcopy(bnds)

    Iend = CartesianIndex(size(dem))

    pits_to_keep = trues(length(pits_))
    # iterate until all interior pits are removed
    # TODO: warn when maxiter is exceeded
    for i=1:maxitr
        n_removed = 0
        for (color, P) in enumerate(pits_)
            copyto!(dir2, dir_) # TODO hack...

            P==CartesianIndex(-1,-1) && continue # already removed pit, skip
            if P.I[1]==1 || P.I[2]==1 || P.I[1]==size(dir,1) || P.I[2]==size(dir,2)
                # don't process pits on the DEM boundary
                continue
            end
            # point on boundary with minimum elevation
            Imin = bnds_[color][findmin(getindex.(Ref(dem), bnds_[color]))[2]]
            @assert c_[Imin]==c_[P]
            # make the outflow and find the next catchment
            min_ = Inf
            target = Imin
            for J in iterate_D9(Imin, dem)
                Imin==J && continue # not the point itself
                c_[J] == color && continue # J is in Imin's catchment
                if dem[J] < min_
                    min_ = dem[J]
                    target = J
                end
            end
            @assert target!=Imin
            # do the flow across the boundary
            _flow_from_to!(Imin, target, dir_, nin_, nout_)

            # reverse directions on path going from Imin to P
            if Imin!=P
                P1 = Imin
                while true
                    # get downstream point
                    P2 = dir2ind(dir2[P1]) + P1 # note usage of `dir2`
                    # The need for dir2 could be removed by first
                    # traversing the flow path without modifying it.
                    # And then do the modification in a second traversal.
                    if P2==P1
                        # With dir2, this should be fixed now.
                        println("This should not happen!")
                        @show color, Imin, P, P1
                        break
                    end
                    # reverse
                    _flow_from_to!(P2, P1, dir_, nin_, nout_)

                    P2==P && break # reached the pit
                    P1 = P2
                end
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
    end

    # remove removed pits
    pits_ = pits_[pits_.!=Ref(CartesianIndex(-1, -1))]
    # redo colors
    cols = Dict([c=>i for (i,c) in enumerate(unique(c_))])
    for i in eachindex(c_)
        c_[i] = cols[c_[i]]
    end

    return dir_, nin_, nout_, pits_, c_, bnds_
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


"Plot DEM, uparea, flow-dir"
function plotit(xs, ys, dem)
    dx2 = step(xs)/2
    area, slen, dir, nout, nin, pits  = waterflows(dem)

    fig, axs = subplots(3,1)
    sca(axs[1])
    contour(xs, ys, Array(dem'))
    colorbar()
    sca(axs[2])
    heatmap(xs, ys, log10.(area))
    sca(axs[3])
    plotdir(xs, ys, dir)
end

function plotdir(xs, ys, dir)
    vecfield = dir2vec.(dir)
    vecfieldx = [v[1] for v in vecfield]
    vecfieldy = [v[2] for v in vecfield]
    quiver(repeat(xs,1, length(ys)), repeat(ys,length(xs),1), vecfieldx, vecfieldy)
end

pits2inds(pits) = ([p.I[1] for p in pits],
                   [p.I[2] for p in pits])
pits2vecs(xs, ys, pits) = (xs[[p.I[1] for p in pits if p!=CartesianIndex(-1,-1)]],
                           ys[[p.I[2] for p in pits if p!=CartesianIndex(-1,-1)]])

"Plot uparea"
function plotarea(xs, ys, dem; prefn=log10)
    area, slen, dir, nout, nin, pits  = waterflows(dem)
    plotarea(xs, ys, area, pits, prefn=prefn)
end
function plotarea(xs, ys, area, pits; prefn=log10)
    px, py = pits2vecs(xs, ys, pits)
    fig, axs = subplots()
    heatmap(xs, ys, prefn.(area))
    scatter(px, py, 1, "w")
end

function plotlakedepth(x, y, dem)
    area, slen, dir, nout, nin, pits  = waterflows(dem)
    demf = fill_dem(dem, pits, dir)
    heatmap(x, y, demf.-dem)
end

function heatmap(xs, ys, mat)
    dx2 = step(xs)/2
    imshow(Array(mat'), origin="lower", extent=(xs[1]-dx2,xs[end]+dx2,ys[1]-dx2,ys[end]+dx2))
    colorbar()
end

end # module
