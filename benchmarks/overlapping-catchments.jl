# julia --project benchmarks/overlapping-catchments.jl
using WhereTheWaterFlows
const WWF = WhereTheWaterFlows
const SG = WWF.Subglacially

function measure(name, f)
    f() # warm up
    GC.gc()
    result = @timed f()
    println(name, ": seconds=", result.time, " allocated_MB=", result.bytes/1e6)
    return result.value
end

n = 256
bed = [1000. - i for i in 1:n, j in 1:n]
surf, source = bed .+ 100., fill(1e-8,n,n)
base = SG.waterflows_subglacial(surf, bed, 200., 0.97, source; gamma=0.)
sinks = base.routing.sinks
groups = [sinks[i:198:end] for i in 1:198]
push!(groups, [CartesianIndex(n÷2,n÷2)]) # interior catchment overlaps terminal groups
out = measure("routing with groups", () -> SG.waterflows_subglacial(surf, bed, 200.,
    0.97, source; gamma=0., ctch_sinks=groups))
masks, a = out.sink_catchments.masks, out.routing.area
old_fluxes = measure("original flux expressions alone", () -> (
    [sum(((source .+ a.dissipation_melt_rate .+ a.pressure_melt_rate)*200^2)[c]) for c in masks],
    [sum((a.dissipation_melt_rate*200^2)[c]) for c in masks],
    [sum((a.pressure_melt_rate*200^2)[c]) for c in masks]))
@assert old_fluxes[1] ≈ out.sink_catchments.fluxes.total
counts = zeros(Float16,n,n,length(groups))
function old_update!(counts, masks)
    for (i,c) in enumerate(masks)
        counts[:,:,i] .+= c
    end
end
function new_update!(counts, masks)
    for (i,c) in enumerate(masks)
        WWF._accumulate_catchment!(view(counts,:,:,i), c)
    end
end
measure("full-layer aggregation", () -> old_update!(counts, masks))
measure("true-cell aggregation", () -> new_update!(counts, masks))
println("Retained mask payload MB: ", sum(sizeof, masks)/1e6)
