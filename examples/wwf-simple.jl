using WhereTheWaterFlows
using Random; Random.seed!(42)

const WWF = WhereTheWaterFlows

n = 120
x = range(-pi, pi, length=n)
dem = sin.(x) .* cos.(x') .+ 0.03*randn(n, n)

out = WWF.waterflows(dem)

println("WhereTheWaterFlows simple example")
println("grid size: ", size(dem))
println("number of sinks: ", length(out.sinks))
println("number of pits after draining: ", length(out.pits))
println("max upslope area: ", maximum(out.area))
