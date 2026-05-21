using WhereTheWaterFlows
using Random

const WWF = WhereTheWaterFlows

Random.seed!(42)

n = 120
xs = range(-pi, pi, length=n)
dem = sin.(xs) .* cos.(xs') .+ 0.03 * randn(n, n)

out = WWF.waterflows(dem)

println("WWF simple example")
println("grid size: ", size(dem))
println("number of sinks: ", length(out.sinks))
println("number of pits after draining: ", length(out.pits))
println("max upslope area: ", maximum(out.area))
