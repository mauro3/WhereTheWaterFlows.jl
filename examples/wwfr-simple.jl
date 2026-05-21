using WhereTheWaterFlows
using Random

const WWFR = WhereTheWaterFlows.Randomly

Random.seed!(42)

n = 100
dx = 100.0
xs = range(-pi, pi, length=n)
dem = sin.(xs) .* cos.(xs')
mask = trues(size(dem))

source = fill(1.0 / dx^2, size(dem))
source_uc = WWFR.Uncertainty(absuc=0.2 / dx^2, reluc=0.15, correlation_length=8.0 * dx)
dem_uc = WWFR.Uncertainty(absuc=0.0, reluc=0.0, correlation_length=dx)

ctch_sinks = [CartesianIndices((2:2, 2:n-1))[:]]

model, sample, reduce! = WWFR.make_fns_subaerial(dx,
                                                 dem, dem_uc,
                                                 source, source_uc,
                                                 ctch_sinks,
                                                 mask)

aggr = WWFR.map_mc(model, sample, reduce!, 8)

println("WWFR simple example (WWF + uncertain source)")
println("samples: ", aggr.n_samples[])
println("mean upslope area (max): ", maximum(aggr.areas_total))
println("mean stream length (max): ", maximum(aggr.stream_length))
println("outlet catchment probability (max): ", maximum(aggr.catchments[:, :, 1]))
