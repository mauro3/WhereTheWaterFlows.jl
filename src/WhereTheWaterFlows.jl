# TODO:
# - check with GRF that it works most of the time
# - better barrier breach algo: select point which has minimal elevation of Imin and target
#   --> what if DEM has a spillway which is a flat area?  Then this method will fail!
# - deinstall packages: Infiltrator, Primes, FFTW

module WhereTheWaterFlows

using StaticArrays, Infiltrator
using OffsetArrays: Origin, no_offset_view

export waterflows

const I11 = CartesianIndex(1,1)
const I22 = CartesianIndex(2,2)
const PIT = 5        # Direction number indicating no flow, this is a "pit", i.e. a
                     # local minimum or a cell in a completely flat area.
const SINK = 10      # A cell where water disappears, typically located at the domain boundary
                     # (when setting bnd_as_sink=true) adjacent to NaN-cells of the DEM.
                     # Note: PITs can also be sinks if drain_pits==false.
const BARRIER = 11   # Direction number indicating no flow into or out of this cell.  All DEM
                     # cells which are NaNs map to this.

"""
Direction numbers.  E.g. dirnums[1,1] will return the number
corresponding to the direction top-left.

Note, I use the conversion that the x-axis corresponds to the row of
the matrix and the y-axis the columns.  To print them in this "normal"
coordinate system use `showme`
"""
const dirnums = SMatrix{3,3}(reverse([ 7   8  9
                                       4  PIT 6
                                       1   2  3]',
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

"""
Translate a D8 direction number into a CartesianIndex (i.e. a flow vector).
Maps dir==BARRIER and dir==SINK to CartesianIndex(0,0) also.
"""
function dir2ind(dir, allow_BARRIER=false)
    allow_BARRIER==false && dir==BARRIER && error("Cannot make CartisianIndex for BARRIER")
    return (dir==SINK || dir==BARRIER) ? CartesianIndex(0,0) : cartesian[dir]
end

"Translate a D8 direction number into a 2D vector."
dir2vec(dir, allow_BARRIER=false) = [dir2ind(dir, allow_BARRIER).I...]

"""
Tests whether a cell `J` with flowdir `dirJ` flows into cell `I`.
"""
flowsinto(J::CartesianIndex, dirJ::Integer, I::CartesianIndex) = ind2dir(I-J) == dirJ

"Return CartesianIndices corresponding to the 8 neighbors and the point itself"
iterate_D9(I, Iend, I1=I11) = max(I1, I-I1):min(Iend, I+I1)
function iterate_D9(I, ar::AbstractMatrix)
    R = CartesianIndices(ar)
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
    d8dir_feature(dem, bnd_as_sink, nan_as_sink)

D8 directions of a DEM and drainage features (nin & nout).

Elevations with NaN map to dir==BARRIER, cells around them will be set to SINK if nan_as_sink==true
or receive no special treatment otherwise.

The argument `bnd_as_sink` determines whether cells at the domain boundary act as sinks.

Return
- dir  - direction, encoded as `dirnums`
- nout - number of outflow cells of a cell (0 or 1)
- nin  - number of inflow cells of a cell (0-8)
- sinks - location of sinks as a `Vector{CartesianIndex{2}}` (sorted) (here dir==SINK)
- pits - location of pits as a `Vector{CartesianIndex{2}}` (sorted) (here dir==PIT)
- dem - DEM, unchanged
- flowdir_extra_output -- nothing (not used by this function, but could be by custom ones)
"""
function d8dir_feature(dem, bnd_as_sink, nan_as_sink)
    # outputs
    dir = fill!(similar(dem, Int8), 0)
    nout = fill!(similar(dem, Bool), false)
    nin = fill!(similar(dem, Int8), 0)
    sinks = CartesianIndex{2}[]
    pits = CartesianIndex{2}[]

    R = CartesianIndices(size(dem))
    Iend = last(R)

    # get dir for all points
    for I in R
        ele = dem[I]

        # always make NaN points BARRIERS
        if isnan(ele)
            dir[I] = BARRIER
            continue
        end
        # make sinks on boundary if bnd_as_sink is set
        if bnd_as_sink && on_outer_boundary(dem,I)
            # make it a sink
            dir[I] = SINK
            continue
        end

        delta_ele = 0.0 # keeps track of biggest elevation change
        dir_ = PIT
        for J in iterate_D9(I, Iend)
            I==J && continue
            ele2 = dem[J]
            if isnan(ele2)
                if nan_as_sink
                    dir_ = SINK
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
        dir[I] = dir_
    end
    # flow features (i.e. nout, nin)
    # (not thread-safe because of push!)
    for I in R
        for J in iterate_D9(I, Iend)
            J==I && continue
            nin[I] += flowsinto(J, dir[J], I)
        end
        if dir[I]==PIT
            nout[I] = false
            push!(pits, I)
        elseif dir[I]==SINK
            nout[I] = false
            push!(sinks, I)
        elseif dir[I]==BARRIER
            nout[I] = false
            # check that nothing flows into it
            @assert nin[I]==0 "Inflow into BARRIER cell detected!  Bug in this code."
        else
            nout[I] = true
        end
    end
    pits_ = Origin(length(sinks)+1)(pits)
    return dir, nout, nin, sinks, pits_, dem, nothing
end


# TODO: think about: if bnd_as_sink routes water along boundary edges first. This would probably
# substantially reduce the number of catchments, as currently every boundary point
# is a pit and thus a catchment (if bnd_as_sink==true).
"""
    waterflows(dem, cellarea=fill!(similar(dem),1), flowdir_fn=d8dir_feature;
               feedback_fn=nothing, drain_pits=true, bnd_as_sink=false)

Does the water flow routing according the D8 algorithm.  Locations of the `dem`
with `NaN`-value are ignored.

args:
- dem -- the DEM (or hydro-potential); array
- cellarea=cellarea=fill!(similar(dem),1) -- the source per cell, defaults to 1.
     - do not use negative values
     - If using physical units then use a volumetric flux per cell, e.g. m3/s.
     - Alternatively, `cellarea` can be a tuple of arrays. Then they are treated/routed
       separately, for instance `(water, tracer)`.  All quantities need to be extensive
       (i.e. additive, e.g. use internal energy and not temperature)
- flowdir_fn=d8dir_feature -- the routing function.  Defaults to the built-in `d8dir_feature`
                              function but could be customized

kwargs:
- feedback_fn -- function which is applied to area-value(s) at each cell once all water
                 of the cell has been accumulated but before the water is routed further downstream.
                 Signature `(uparea, ij, dir) -> new_uparea`
- drain_pits -- whether to route through pits (true)
- bnd_as_sink (true) -- whether the domain boundary should be sinks, i.e. adjacent cells
                 can drain into them, or whether to ignore them.
- nan_as_sink (true) -- whether NaN cells in the DEM should make adjacent cells a sink.
- stacksize (2^13 * 2^10) -- size of the call-stack in _flowrouting_catchments!, which is prone to
                 StackOverflowError.  Note however, that OutOfMemory errors are likely if increased.


Returns
- area -- upslope area (or a tuple of upslope areas if cellarea is a tuple too)
- stream-length -- length of stream to the farthest source (number of cells traversed)
- dir -- flow direction at each location
- nout -- whether the point has outlflow.  I.e. nout[I]==0 --> I is a pit
- nin -- number of inflow cells
- sinks -- location of sinks as Vector{CartesianIndex{2}}
- pits -- location of pits as Vector{CartesianIndex{2}}
- c -- catchment map (color numbers ∈ 1:length(sinks) are for sinks, others for pits)
- bnds -- boundaries between catchments.  The boundary to the exterior/NaNs is not in here.
- flowdir_extra_output -- extra output of the flowdir_fn, which is `nothing` for the default
"""
function waterflows(dem, cellarea=fill!(similar(dem),1), flowdir_fn=d8dir_feature;
                    feedback_fn=nothing, drain_pits=true, bnd_as_sink=true, nan_as_sink=true, stacksize=2^13 * 2^10)
    if drain_pits && !bnd_as_sink
        if !nan_as_sink || (nan_as_sink && sum(isnan.(dem))==0)
            error("No sinks in the domain.  Consider setting bnd_as_sink and/or nan_as_sink and add NaNs to the `dem`.")
        end
    end

    dir, nout, nin, sinks, pits, dem4drainpits, flowdir_extra_output = flowdir_fn(dem, bnd_as_sink, nan_as_sink)
    area, slen, c = flowrouting_catchments(dir, sinks, pits, cellarea, feedback_fn, stacksize)
    bnds = make_boundaries(c, collect(axes(pits)[1]))
    if drain_pits
        drainpits!(dir, nin, nout, sinks, pits, c, bnds, dem4drainpits)
        area, slen, c = flowrouting_catchments(dir, sinks, pits, cellarea, feedback_fn, stacksize)
    end
    #area[isnan.(dem)] .= NaN
    return area, slen, dir, nout, nin, sinks, pits, c, bnds, flowdir_extra_output
end

"""
    flowrouting_catchments(dir, pits, cellarea, feedback_fn, stacksize)

Recursively calculate flow-routing and catchments from
- dir - direction field
- pits - pit coordinates
- cellarea - water input
- stacksize - how large the call-stack is (in bytes)

Return:
- upslope area
- stream length (number of cells traversed)
- catchments Matrix{Int}.  Value==0 corresponds to NaNs in the DEM
  which are not pits (i.e. where no water flows into).

Note: this function may cause a stackoverflow on very big catchments.
"""
function flowrouting_catchments(dir, sinks, pits, cellarea, feedback_fn, stacksize) # on linux standard is 2^13 * 2^10
    c = fill!(similar(dir, Int), 0) # catchment color of BARRIER is 0
    slen = fill!(similar(dir, Int), 0)
 @show   np = length(pits)
    # Some setup to allow both array and tuple-of-array inputs for cellarea and feedback_fn:
    area = init_area(dir, cellarea) # (matrix,) or tuple of matrices
    cellarea_ = cellarea isa Tuple ? cellarea : (cellarea, )
    feedback_fn_ = if feedback_fn===nothing
        nothing
    else
        cellarea isa Tuple ? feedback_fn : (uparea, ij, dir) -> (feedback_fn(uparea[1], ij, dir),)
    end

    # recursively traverse the drainage tree in up-flow direction,
    # starting at all sinks and pits (the pits will likely be removed with drainpits! eventually)
    #Threads.@threads
    for color = 1:length(sinks)
        sink = sinks[color]

        # This is a dirty trick to increase the call-stack size
        # https://stackoverflow.com/questions/71956946/how-to-increase-stack-size-for-julia-in-windows
        # Note that stacksize is a undocumented argument of Task!
        wait(schedule( Task(() -> _flowrouting_catchments!(area, slen, c, dir, cellarea_, feedback_fn_, color, sink),
            stacksize) ))
    end
    for color in axes(pits)[1]
       pit = pits[color]

        # This is a dirty trick to increase the call-stack size
        # https://stackoverflow.com/questions/71956946/how-to-increase-stack-size-for-julia-in-windows
        # Note that stacksize is a undocumented argument of Task!
        wait(schedule( Task(() -> _flowrouting_catchments!(area, slen, c, dir, cellarea_, feedback_fn_, color, pit),
                            stacksize) ))
    end
    # unwrap tuple if cellarea was not a tuple:
    if !(cellarea isa Tuple)
        area = area[1]
    end
    return area, slen, c
end
# initialize output cellarea:
# - as tuple or one array
# - as tuple, if cellarea is a tuple
init_area(dir, cellarea) = (fill!(similar(dir,Float64), 0), )
init_area(dir, cellarea::Tuple) = map(x -> fill!(similar(dir,Float64), 0), cellarea)

# modifies c and area
function _flowrouting_catchments!(area, len, c, dir, cellarea, feedback_fn, color, ij)
    # assign catchment (solely dependent on `dir`)
    c[ij] = color
    n = length(cellarea)

    # proc upstream points
    slen = 0 # note: solely dependent on `dir`
    uparea = getindex.(cellarea, Ref(ij))
    slen = max(slen, 1)
    for IJ in iterate_D9(ij, c)
        if ij==IJ
            continue
        elseif flowsinto(IJ, dir[IJ], ij)
            uparea_, slen_ = _flowrouting_catchments!(area, len, c, dir, cellarea, feedback_fn, color, IJ)
            uparea = uparea .+ uparea_
            slen = max(slen, slen_+1) # TODO take diagonal into account
        end
    end
    # feedback of areas with each other and onto themselves
    if feedback_fn!==nothing
        uparea = feedback_fn(uparea, ij, dir)
    end
    setindex!.(area, uparea, Ref(ij)) # the setindex! is needed for the broadcasting over the area-tuple to work
    len[ij] = slen
    return uparea, slen
end

"""
    make_boundaries(catchments, pit_colors)

Make vectors of boundary cells for catchments of `pit_colors` (as the name suggests,
typically just the pit-catchment colors.)  Assumes that pit_colors is sorted.

Note that cells along the edge of the domain[1] as well as cells only bordering
BARRIER-cells are not included (because the algorithms do not need to traverse them).

Return:
- bnds -- Vector of Vector{CartesianIndex{2}} containing the cells which are
          on the boundary of said catchment.

[1] Note that when bnd_as_sink==true then no cells along the boundary will belong
    to a pit-catchment.
"""
function make_boundaries(catchments, pit_colors)
    length(pit_colors)==0 && return Vector{CartesianIndex{2}}[]
    pc1 = pit_colors[1]
    bnds = Origin(pc1)([CartesianIndex{2}[] for i in 1:length(pit_colors)])
    for R in CartesianIndices(axes(catchments))
        c = catchments[R]
        c < pc1 && continue # don't find boundaries for sink catchments
        bnd = bnds[c]
        for I in iterate_D9(R, catchments)
            I==R && continue
            co = catchments[I]
            if co!=c && co!=0 # co!=0 means its a BARRIER cell
                push!(bnd, R)
                break
            end
        end
    end
    return bnds
end

"""
    _prune_boundary!(del, bnds, catchments::AbstractMatrix, color, colormap)

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

function drainpits_old!(dir, nin, nout, sinks, pits, c, bnds, dem)
    length(pits)==0 && return nothing
    nsinks = length(sinks)
    ncolors = length(sinks) + length(pits)

    maxiter = 100
    Iend = CartesianIndex(size(dem))

    pits_to_delete = Origin(nsinks+1)(falses(length(pits)))

    # Table which translates the old pit-color to the new pit-color.
    # Note that only the currently processed color might change.
    # Initialize to the id-map (0 is used for off-points)
    colormap = collect(1:ncolors)

    #colormap_inv = [[i] for i=1:ncolors]

    # work array as otherwise much time is spent allocating
    # (if below loop is ever multi-threaded, this will need one per thread)
    del_workarray = Int[]
    sizehint!(del_workarray, maximum(length.(bnds)))

    # iterate until all interior pits are removed (not really needed, I think)
    @inbounds for iter=1:maxiter
        n_removed = 0
        for pit_color in axes(pits)[1]
            P = pits[pit_color]
            # Already removed pit, skip
            P==CartesianIndex(-1,-1) && continue

            @assert colormap[pit_color] == pit_color "$(colormap[pit_color]) $(pit_color)"

            # # If there are no more boundaries left, stop.  This should only occur when
            # # bnd_as_sink==false and when there are only interior pits.
            # if all((isempty(b) for b in bnds))
            #     no_drainage_across_boundary = true # set flag to warn later
            #     break
            # end

            # If there are no boundaries for this point error
            isempty(bnds[pit_color]) && error("Found a pit-catchment which has no outflow.  If this is intended consider setting the DEM-cell to a NaN at $P and setting nan_as_sink.")

            # Find point on catchment boundary with minimum elevation
            minn = typemax(eltype(dem))
            Imin = CartesianIndex(-1,-1)
            for I in bnds[pit_color]
                if dem[I] < minn
                    minn = dem[I]
                    Imin = I
                    dir[I]==BARRIER && error("Bug in code")
                end
            end

            @assert colormap[c[Imin]]==c[P]==pit_color

            # Make the outflow and find the next catchment
            min_ = Inf
            target = Imin # target is the cell in the next catchment
            for J in iterate_D9(Imin, dem)
                Imin==J && continue # do not look at the point itself
                cc = c[J]
                cc = cc==0 ? 0 : colormap[cc] # to catch BARRIER cells
                cc == pit_color && continue # J is in Imin's catchment
                if dem[J] < min_
                    min_ = dem[J]
                    target = J
                end
            end
            if target==Imin
                error("Bug in Program")
            end

            # Reverse directions on path going from Imin to P
            P1 = Imin
            P2 = Imin + dir2ind(dir[Imin]) # do lookup ahead of:
            _flow_from_to!(Imin, target, dir, nin, nout) # first, do the flow across the boundary
            while P1!=P
                # get down-downstream point (do ahead lookup as a undisturbed dir is needed)
                P3 = P2 + dir2ind(dir[P2])
                # reverse
                _flow_from_to!(P2, P1, dir, nin, nout)

                P1 = P2
                P2 = P3
            end

            # remove from list of pits
            pits_to_delete[pit_color] = true
            pits[pit_color] = CartesianIndex(-1,-1)

            # update colormap and boundaries
            othercolor = colormap[c[target]]
            ## this is slow
            for i in eachindex(colormap)
                if pit_color==colormap[i]
                    colormap[i] = othercolor
                end
            end
            # ## thus do this instead:
            # if c[target]>nsinks
            #     append!(colormap_inv[othercolor], colormap_inv[pit_color])
            #     empty!(colormap_inv[pit_color])
            #     for col in colormap_inv[othercolor]
            #         colormap[col] = othercolor
            #     end
            # end

            if othercolor>nsinks
                append!(bnds[othercolor], bnds[pit_color])
            end
            empty!(bnds[pit_color])
            if othercolor>nsinks
                _prune_boundary!(del_workarray, bnds, c, othercolor, colormap)
            end
            n_removed +=1
        end
        n_removed==0 && break
        if iter==maxiter
            error("Maximum number of iterations reached in `drainpits`: $i")
        else
            iter>1 && @show iter
        end
    end

    # # make new colors consecutive
    # newcolors = unique(colormap)
    # d = Dict((c=>i for (i,c) in enumerate(newcolors)))
    # colormap2 = [d[c] for c in colormap]

    # remove removed pits and bnds
    deleteat!(no_offset_view(pits), no_offset_view(pits_to_delete))
    deleteat!(no_offset_view(bnds), no_offset_view(pits_to_delete))

    # # update catchments -> done with `flowrouting_catchments`
    return nothing
end

"""
    drainpits!(dir, nin, nout, sinks, pits, c, bnds, dem)

Update in place the direction field such that it drains pits. This is done by
reversing the flow connecting the lowest point (which can drain) on the catchment
boundary to the pit for each such catchment.

Update in place dir, nin, nout, pits (sorted), c, bnds

TODO: this is the performance bottleneck.
"""
function drainpits!(dir, nin, nout, sinks, pits, c, bnds, dem)
    length(pits)==0 && return nothing
    nsinks = length(sinks)
    ncolors = length(sinks) + length(pits)

    maxiter = 100
    Iend = CartesianIndex(size(dem))

    pits_to_delete = Origin(nsinks+1)(falses(length(pits)))

    # Table which translates the old pit-color to the new pit-color.
    # Note that only the currently processed color might change.
    # Initialize to the id-map (0 is used for off-points)
    colormap = collect(1:ncolors)

    colormap_inv = [[i] for i=1:ncolors]

    # work array as otherwise much time is spent allocating
    # (if below loop is ever multi-threaded, this will need one per thread)
    del_workarray = Int[]
    sizehint!(del_workarray, maximum(length.(bnds)))

    # iterate until all interior pits are removed (not really needed, I think)
    @inbounds for iter=1:maxiter
        n_removed = 0
        for pit_color in axes(pits)[1]
            P = pits[pit_color]
            # Already removed pit, skip
            P==CartesianIndex(-1,-1) && continue

            @assert colormap[pit_color] == pit_color "$(colormap[pit_color]) $(pit_color)"

            # # If there are no more boundaries left, stop.  This should only occur when
            # # bnd_as_sink==false and when there are only interior pits.
            # if all((isempty(b) for b in bnds))
            #     no_drainage_across_boundary = true # set flag to warn later
            #     break
            # end

            # error if there are no boundaries for this cell
            isempty(bnds[pit_color]) && error("Found a pit-catchment which has no outflow.  If this is intended consider setting the DEM-cell to a NaN at $P and setting nan_as_sink.")

            # Find cell on catchment boundary with minimum spillway elevation
            spillway = typemax(eltype(dem))  # elevation of spillway
            P_i = CartesianIndex(-1,-1) # cell on spillway within the pit_color-catchment
            P_o = CartesianIndex(-1,-1) # outside cell of spillway
            for I in bnds[pit_color]
                spillway_I = dem[I] # spillway at cell I is at this elevation or higher
                spillway_I > spillway && continue
                # check cells on other catchment(s) to see if spillway_I has to be increased:
                J_min = CartesianIndex(-1,-1) # lowest cell in other catchment
                ele_min = typemax(eltype(dem))# elevation of that cell
                for J in iterate_D9(I, dem)
                    I==J && continue # do not look at the cell itself
                    c[J] == 0 && continue # J is a BARRIER cell
                    colormap[c[J]] == pit_color && continue # J is in I's catchment
                    if dem[J] < ele_min
                        ele_min = dem[J]
                        J_min = J
                    end
                end
                J_min == CartesianIndex(-1,-1) && error("Bug in program")
                spillway_I = max(spillway_I, ele_min)
                # check if spillway_I is a new lowest spillway
                if spillway_I < spillway
                    P_i = I
                    P_o = J_min
                    spillway = spillway_I
                end
            end
            @assert pit_color != colormap[c[P_o]] "Error: both spillway cells in the same catchment"

            # Reverse directions on path going from P_o via P_i to P
            P1 = P_i
            P2 = P_i + dir2ind(dir[P_i]) # do lookup ahead of:
            _flow_from_to!(P_i, P_o, dir, nin, nout) # first, do the flow across the boundary
            while P1!=P
                # get down-downstream cell (do ahead lookup as a undisturbed dir is needed)
                P3 = P2 + dir2ind(dir[P2])
                # reverse
                _flow_from_to!(P2, P1, dir, nin, nout)

                P1 = P2
                P2 = P3
            end

            # remove from list of pits
            pits_to_delete[pit_color] = true
            pits[pit_color] = CartesianIndex(-1,-1)

            # update colormap and boundaries
            # (now this catchment becomes the color of the catchment it now overspills into)
            othercolor = colormap[c[P_o]]
            # ## this is slow
            # for i in eachindex(colormap)
            #     if pit_color==colormap[i]
            #         colormap[i] = othercolor
            #     end
            # end
            ## thus do this instead:
            if c[P_o]>nsinks
                append!(colormap_inv[othercolor], colormap_inv[pit_color])
                empty!(colormap_inv[pit_color])
                for col in colormap_inv[othercolor]
                    colormap[col] = othercolor
                end
            end

            if othercolor>nsinks
                append!(bnds[othercolor], bnds[pit_color])
            end
            empty!(bnds[pit_color])
            if othercolor>nsinks
                _prune_boundary!(del_workarray, bnds, c, othercolor, colormap)
            end
            n_removed +=1
        end
        n_removed==0 && break
        if iter==maxiter
            error("Maximum number of iterations reached in `drainpits`: $i")
        else
            iter>1 && @show iter
        end
    end

    # # make new colors consecutive
    # newcolors = unique(colormap)
    # d = Dict((c=>i for (i,c) in enumerate(newcolors)))
    # #colormap2 = [d[c] for c in colormap]

    # remove removed pits and bnds
    deleteat!(no_offset_view(pits), no_offset_view(pits_to_delete))
    deleteat!(no_offset_view(bnds), no_offset_view(pits_to_delete))

    # # update catchments -> done with `flowrouting_catchments`
    return nothing
end

function drainpits2!(dir, nin, nout, sinks, pits, c, bnds, dem)
    length(pits)==0 && return nothing
    @assert length(pits)==length(bnds)

    nsinks = length(sinks)
    npits = length(pits)
    ncolors = nsinks + npits

    tmp = eltype(pits)[]
    sizehint!(tmp, npits÷4)
    new_pits = Origin(nsinks+1)(tmp) # when two re-directed branches meet at a spillway,
                             # they create a new pit at the spillway


    maxiter = 100
    Iend = CartesianIndex(size(dem))

    dirty_catchment = falses(nsinks+npits)
    # pits_to_delete = Origin(nsinks+1)(falses(length(pits)))

    for pit_color in axes(pits)[1]
        P = pits[pit_color]

        # If there are no boundaries for this point error
        isempty(bnds[pit_color]) && error("Found a pit-catchment which has no outflow.  If this is intended consider setting the DEM-cell to a NaN at $P and setting nan_as_sink.")

        # Find cell on catchment boundary with minimum spillway elevation
        spillway = typemax(eltype(dem))  # elevation of spillway
        P_i = CartesianIndex(-1,-1) # cell on spillway within the pit_color-catchment
        P_o = CartesianIndex(-1,-1) # outside cell of spillway
        for I in bnds[pit_color]
            spillway_I = dem[I] # spillway at cell I is at this elevation or higher
            spillway_I > spillway && continue
            # check cells on other catchment(s) to see if spillway_I has to be increased:
            J_min = CartesianIndex(-1,-1) # lowest cell in other catchment
            ele_min = typemax(eltype(dem))# elevation of that cell
            for J in iterate_D9(I, dem)
                I==J && continue # do not look at the cell itself
                c[J] == 0 && continue # J is a BARRIER cell
                c[J] == pit_color && continue # J is in I's catchment
                if dem[J] < ele_min
                    ele_min = dem[J]
                    J_min = J
                end
            end
            J_min == CartesianIndex(-1,-1) && error("Bug in program")
            spillway_I = max(spillway_I, ele_min)
            # check if spillway_I is a new lowest spillway
            if spillway_I < spillway
                P_i = I
                P_o = J_min
                spillway = spillway_I
            end
        end
        # potentially adjust P_i such that it is lowestmost cell bordering on P_o
        ele_min = typemax(eltype(dem))
        for J in iterate_D9(P_o, dem)
            P_o==J && continue # do not look at the cell itself
            c[J] == 0 && continue # J is a BARRIER cell
            c[J]!=c[P_i] && continue # J is not in P_i's catchment
            if dem[J] < ele_min
                ele_min = dem[J]
                P_i = J
            end
        end

        if c[P_o]>nsinks
            # if the other catchment is in a pit-catchment, mark it dirty
            dirty_catchment[c[P_o]] = true
        end

        # if P∈[CartesianIndex(9,37), CartesianIndex(14,37), CartesianIndex(8,33)]
        #     @infiltrate
        # end

        # Reverse directions on path going from P_i to P and
        # create flow from P_i to P_o (i.e. across the spillway)
        P1 = P_i
        P2 = P_i + dir2ind(dir[P_i]) # do lookup ahead of re-direction
        if flowsinto(P_o, dir[P_o], P_i)
            # This means this spillway was already used from the other side.
            # This means a new pit is created at the spillway.
            @assert dirty_catchment[c[P_o]]
            push!(new_pits, P_o)
            allow_P_o_pit = true
        else
            allow_P_o_pit = false
        end
        _flow_from_to!(P_i, P_o, dir, nin, nout, allow_P_o_pit) # first, do the flow across the boundary
        while P1!=P
            # get down-downstream point (do ahead lookup as a undisturbed dir is needed)
            P_next = P2 + dir2ind(dir[P2])
            # reverse
            _flow_from_to!(P2, P1, dir, nin, nout)

            P1,P2 = P2,P_next
        end
        #pits_to_delete[pit_color] = true
    end
    # update pits
    #deleteat!(no_offset_view(pits), no_offset_view(pits_to_delete))

    # bnds is now invalid
    empty!(bnds)
    return new_pits
end

function drainpits2_old!(dir, nin, nout, sinks, pits, c, bnds, dem)
    length(pits)==0 && return nothing
    @assert length(pits)==length(bnds)

    nsinks = length(sinks)
    npits = length(pits)
    ncolors = nsinks + npits

    maxiter = 100
    Iend = CartesianIndex(size(dem))

    dirty_catchment = falses(nsinks+npits)
    pits_to_delete = Origin(nsinks+1)(falses(length(pits)))

    for pit_color in axes(pits)[1]
        P = pits[pit_color]
        # dirty catchments cannot be processed in this round, but need to wait for the next iteration
        dirty_catchment[pit_color] && continue

        # If there are no boundaries for this point error
        isempty(bnds[pit_color]) && error("Found a pit-catchment which has no outflow.  If this is intended consider setting the DEM-cell to a NaN at $P and setting nan_as_sink.")

        # Find point on catchment boundary with minimum elevation
        minn = typemax(eltype(dem))
        Imin = CartesianIndex(-1,-1)
        for I in bnds[pit_color]
            if dem[I] < minn
                minn = dem[I]
                Imin = I
                dir[I]==BARRIER && error("Bug in code")
            end
        end

        # Make the outflow and find the next catchment
        min_ = typemax(eltype(dem))
        target = Imin # target is the cell in the next catchment
        for J in iterate_D9(Imin, dem)
            Imin==J && continue # do not look at the point itself
            cc = c[J]
            cc==0 && continue # do not look at BARRIER cells
            cc == pit_color && continue # J is in Imin's catchment
            if dem[J] < min_
                min_ = dem[J]
                target = J
            end
        end
        target==Imin && error("Bug in Program")
        if c[target]>nsinks
            # if the target is in a pit-catchment, mark that pit-catchment as dirty
            dirty_catchment[c[target]] = true
        end

        # Reverse directions on path going from Imin to P
        P1 = Imin
        P2 = Imin + dir2ind(dir[Imin]) # do lookup ahead of:
        _flow_from_to!(Imin, target, dir, nin, nout) # first, do the flow across the boundary
        while P1!=P
            # get down-downstream point (do ahead lookup as a undisturbed dir is needed)
            P_next = P2 + dir2ind(dir[P2])
            # reverse
            _flow_from_to!(P2, P1, dir, nin, nout)

            P1,P2 = P2,P_next
        end
        pits_to_delete[pit_color] = true
    end
    # update pits
    deleteat!(no_offset_view(pits), no_offset_view(pits_to_delete))

    # bnds is now invalid
    empty!(bnds)
    return nothing
end

"""
Update dir, nin, and nout such that flow at P1 is now from P1 to P2.

It can potentially modify `dir, nin, nout` at three locations
- P1: dir, nout, nin
- P2: nin
  - if allow_P2_pit==true, then P2's dir, nout can also be modified.
- P3 (previous receiver cell of P1): nin

Note that
- if P1 and P2 lie in the same catchment, then P3 is also in that catchment.
- if flow was from P2 to P1, then P2 has to become a pit (PIT) to keep dir consistent
"""
function _flow_from_to!(P1, P2, dir, nin, nout, allow_P2_pit=false)
    # already right
    ind2dir(P2-P1)==dir[P1] && return nothing
    # if one is a BARRIER the error
    (dir[P1]==BARRIER || dir[P2]==BARRIER) && error("Bug: trying to re-route into/out-of BARRIER point.")
    dir[P1]==SINK && error("Bug: trying to re-route flow out-of SINK point.")

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
        if allow_P2_pit
            dir[P2] = PIT
            nout[P2] = false
            nin[P1] -= 1
        else
            error("P2 $P2 would need to become a pit for P1 $P1 to flow into it (if desired, set option `allow_P2_pit`).")
        end
    end
    return nothing
end


## Post processing
include("postproc.jl")

## Plotting via Plt
struct Plt end

"""
Plotting functions can be accessed via `plt` after loading PyPlot.

Example
```
WhereTheWaterFlows.plt.plotit(dem)
```
"""
plt = Plt()
function Base.getproperty(::Plt, name::Symbol)
    ext = Base.get_extension(@__MODULE__, :PyPlotExt)
    if isnothing(ext)
        error("Need to `using PyPlot` to make plotting available")
    end
    return getproperty(ext, name)
end

function Base.propertynames(::Plt)
    ext = Base.get_extension(@__MODULE__, :PyPlotExt)
    if isnothing(ext)
        error("Need to `using PyPlot` to make plotting available")
    end
    return propertynames(ext)
end

end # module
