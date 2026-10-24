# Run in an environment containing BenchmarkTools.
using BenchmarkTools, Random

rng = Xoshiro(42)
mask = falses(10_000, 1_000)
# Construct exact byte-level clustering; no unused tail bits for this shape.
bytes = reinterpret(UInt8, mask.chunks)
order = randperm(rng, length(bytes))
nfull = 900_000 ÷ 8
bytes[order[1:nfull]] .= 0xff
for j in nfull+1:nfull+100_000
    bytes[order[j]] = UInt8(1) << rand(rng, 0:7)
end
@assert count(mask) == 1_000_000
@assert count(==(0xff), bytes) == nfull
layer = zeros(Float16, size(mask))

function accumulate_normal!(layer, mask)
    for i in eachindex(layer, mask)
        if mask[i]
            layer[i] += one(eltype(layer))
        end
    end
    return layer
end

function accumulate_findnext!(layer, mask)
    bits = vec(mask)
    i = findnext(bits, 1)
    while !isnothing(i)
        layer[i] += one(eltype(layer))
        i = findnext(bits, i + 1)
    end
    return layer
end

function accumulate_logical!(layer, mask)
    for i in Base.LogicalIndex{Int}(vec(mask))
        layer[i] += one(eltype(layer))
    end
    return layer
end

println("Julia ", VERSION, "; CPU ", Sys.CPU_NAME, "; threads ", Threads.nthreads())
println("10000×1000 BitMatrix; 10% true, 90% of true bits in full bytes; Float16; seed 42")
for f in (accumulate_normal!, accumulate_findnext!, accumulate_logical!)
    fill!(layer, 0)
    @assert f(layer, mask) === layer
    @assert layer == mask
    println(f)
    # Reset outside timing; one update per sample avoids accumulating past 2048.
    @btime $f($layer, $mask) setup=(fill!($layer, 0)) evals=1 seconds=2
    @assert layer == mask
end
