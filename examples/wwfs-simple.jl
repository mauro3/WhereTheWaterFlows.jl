using WhereTheWaterFlows
using Random

const WWFS = WhereTheWaterFlows.Subglacially

Random.seed!(42)

n = 140
dx = 100.0
x = range(0, step=dx, length=n)
y = x

X = repeat(collect(x), 1, n)
Y = repeat(collect(y)', n, 1)

surfdem = 1200.0 .+ 0.015 .* X .+ 0.015 .* Y .+ 8.0 .* randn(n, n)
beddem = 950.0 .+ 0.01 .* X .+ 0.01 .* Y .- 120.0 .* exp.(-((X .- x[end] / 2).^2 .+ (Y .- y[end] / 2).^2) ./ (2 * 5000.0^2))
surfdem = max.(surfdem, beddem .+ 10.0)

out = WWFS.waterflows_subglacial(surfdem, beddem, dx; gamma=WWFS.GAMMA)

println("WWFS simple example")
println("grid size: ", size(surfdem))
println("number of sinks: ", length(out.routing.sinks))
println("supercooling-cell count: ", sum(out.pressmelt.sc_locs))
println("max subglacial discharge proxy: ", maximum(out.routing.area.total))
