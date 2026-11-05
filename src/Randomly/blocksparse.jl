# AI disclosure: developed with an AI coding assistant; requires maintainer review.

using BlockArrays: BlockArray, blocks, mortar
using FillArrays: Zeros

const _CATCHMENT_BLOCKSIZE = 50

function _catchment_blocksize(blocksize)
    (blocksize isa Integer && blocksize > 0) || throw(ArgumentError(
        "catchment_blocksize must be a positive integer"))
    return Int(blocksize)
end

"""
    _init_catchments(sz, ngroups, storage; blocksize=_CATCHMENT_BLOCKSIZE)

Create Float16 catchment counts with `:blocksparse` or `:dense` storage.

Sparse storage uses square `blocksize × blocksize × 1` blocks, with smaller
spatial edge blocks. Empty blocks share lazy zeros; each occupied block is
allocated independently on first use. `blocksize` defaults to 50 and is ignored
for `:dense` storage.
"""
function _init_catchments(sz, ngroups, storage; blocksize=_CATCHMENT_BLOCKSIZE)
    storage === :dense && return zeros(Float16, sz..., ngroups)
    storage === :blocksparse || throw(ArgumentError(
        "catchment_storage must be :blocksparse or :dense"))
    blocksize = _catchment_blocksize(blocksize)
    rows = [min(blocksize, sz[1]-i+1) for i in 1:blocksize:sz[1]]
    cols = [min(blocksize, sz[2]-j+1) for j in 1:blocksize:sz[2]]
    tiles = Array{AbstractArray{Float16,3}}(undef, length(rows), length(cols), ngroups)
    emptytiles = Dict{Tuple{Int,Int},AbstractArray{Float16,3}}()
    for I in CartesianIndices(tiles)
        nr, nc = rows[I[1]], cols[I[2]]
        tiles[I] = get!(emptytiles, (nr, nc)) do
            Zeros{Float16}(nr, nc, 1)
        end
    end
    return mortar(tiles, rows, cols, ones(Int, ngroups))
end

"""
    _accumulate_catchment!(layer, mask)
    _accumulate_catchment!(counts, mask, group)

Add one at each true mask cell, allocating only newly occupied sparse blocks.
"""
function _accumulate_catchment!(layer, mask)
    size(layer) == size(mask) || throw(DimensionMismatch("Aggregate and mask sizes differ"))
    for i in Base.LogicalIndex{Int}(vec(mask))
        layer[i] += one(eltype(layer))
    end
    return layer
end

function _accumulate_catchment!(counts::AbstractArray{T,3}, mask, group) where T
    _accumulate_catchment!(view(counts, :, :, group), mask)
    return counts
end

function _accumulate_catchment!(counts::BlockArray{T,3}, mask, group) where T
    size(counts)[1:2] == size(mask) || throw(DimensionMismatch("Aggregate and mask sizes differ"))
    checkbounds(counts, :, :, group)
    tiles = blocks(counts)
    nx = size(counts, 1)
    # Derive the regular block size from the actual block layout. The first
    # tile is full-sized unless the whole domain is smaller than one block, in
    # which case there is only a single block in that dimension anyway.
    firsttile = tiles[1, 1, group]
    blockrows, blockcols = size(firsttile, 1), size(firsttile, 2)
    previous = (-1, -1)
    tile = nothing
    # Keep the packed-bit path for Raster/DimArray-wrapped masks as well.
    for k in Base.LogicalIndex{Int}(vec(mask))
        col, row = divrem(k-1, nx)
        bi, ii = divrem(row, blockrows)
        bj, jj = divrem(col, blockcols)
        if (bi, bj) != previous
            tile = tiles[bi+1, bj+1, group]
            if tile isa Zeros
                tile = zeros(T, size(tile))
                tiles[bi+1, bj+1, group] = tile
            end
            previous = (bi, bj)
        end
        @inbounds (tile::Array{T,3})[ii+1, jj+1, 1] += one(T)
    end
    return counts
end

"""
    _normalize_catchments!(counts, n)

Divide catchment counts in place, leaving lazy zero blocks unallocated.
"""
_normalize_catchments!(counts, n) = (counts ./= n)
function _normalize_catchments!(counts::BlockArray{T,3}, n) where T
    for tile in blocks(counts)
        tile isa Zeros || _normalize_catchment_tile!(tile, n)
    end
    return counts
end

function _normalize_catchment_tile!(tile::Array{T,3}, n) where T
    tile ./= n
    return nothing
end
