# AI disclosure: developed with an AI coding assistant; requires maintainer review.
# Run: julia --project=. benchmarks/blocksparse-catchments.jl
using WhereTheWaterFlows, Statistics
const WWF = WhereTheWaterFlows
const WWFS = WWF.Subglacially
const WWFR = WWF.Randomly

function update!(counts,masks)
    for (group,mask) in enumerate(masks)
        WWFR._accumulate_catchment!(counts,mask,group)
    end
    return counts
end

function compare_storage(width)
    nx,ny,ngroups=1024,1024,24
    masks=[falses(nx,ny) for _ in 1:ngroups]
    for (g,mask) in enumerate(masks)
        i,j=1+mod(137g,nx-width),1+mod(97g,ny-width)
        mask[i:i+width-1,j:j+width-1] .= true
    end
    reference=nothing
    println("\nSynthetic $(width)×$(width) masks: $(nx)×$(ny), $ngroups overlapping groups")
    for storage in (:dense,:blocksparse)
        small=WWFR._init_catchments((3,3),1,storage)
        WWFR._accumulate_catchment!(small,trues(3,3),1)
        WWFR._normalize_catchments!(small,8)
        setup=@timed WWFR._init_catchments((nx,ny),ngroups,storage)
        counts=setup.value
        update!(counts,masks) # warm update kernel; first of eight draws
        trials=[@timed update!(counts,masks) for _ in 1:7]
        normalize=@timed WWFR._normalize_catchments!(counts,8)
        if reference===nothing
            reference=Array(counts)
        else
            @assert Array(counts)==reference
        end
        println(storage,": object MiB=",round(Base.summarysize(counts)/1024^2;digits=3),
            "; init ms=",round(setup.time*1000;digits=3),
            "; median update ms=",round(median(t.time for t in trials)*1000;digits=3),
            "; median update bytes=",median(t.bytes for t in trials),
            "; normalize ms=",round(normalize.time*1000;digits=3),
            "; normalize bytes=",normalize.bytes)
    end
end

compare_storage(128)
GC.gc()
compare_storage(512)
