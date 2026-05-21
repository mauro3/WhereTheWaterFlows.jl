using WhereTheWaterFlows
using Random; Random.seed!(42)

const WWFS = WhereTheWaterFlows.Subglacially

# artifical topography
n = 140
dx = 100.0
y = x = range(0, step=dx, length=n)
surfdem = @. 1200 + 0.015*x + 0.015*y' + 8*$randn(n, n)
beddem = @. 950 + 0.01*x + 0.01*y' - 120*exp(-((x - x[end]/2)^2 + (y' - y[end]/2)^2) / (2*5000^2))
surfdem = max.(surfdem, beddem .+ 10.0)

out = WWFS.waterflows_subglacial(surfdem, beddem, dx; gamma=WWFS.GAMMA)

println("WhereTheWaterFlowsSubglacially simple example (WWF + Shreve-potential based routing)")
println("grid size: ", size(surfdem))
println("number of sinks: ", length(out.routing.sinks))
println("supercooling-cell count: ", sum(out.pressmelt.sc_locs))
println("max subglacial discharge proxy: ", maximum(out.routing.area.total))
