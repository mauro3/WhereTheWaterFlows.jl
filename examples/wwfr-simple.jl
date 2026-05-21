using WhereTheWaterFlows
using Random, Statistics; Random.seed!(42)

const WWFR = WhereTheWaterFlows.Randomly

n = 100
dx = 100.0
x = range(-pi, pi, length=n)
dem = sin.(x) .* cos.(x') .+ 0.03*randn(n, n)
source = fill(1e-3/dx^2, size(dem))

# Uncertainties from which the Monte Carlo sampler will draw
source_uc = WWFR.Uncertainty(absuc=0, reluc=0.25, correlation_length=n/5*dx) # uncertainty in source
dem_uc = WWFR.Uncertainty() # no uncertainty on DEM
ctch_sinks = [CartesianIndices((2:2, 2:n-1))[:]] # aggregate all discharge going into the left domain margin

model, sample, reduce! = WWFR.make_fns_subaerial(dx,
                                                 dem, dem_uc,
                                                 source, source_uc,
                                                 ctch_sinks)

aggr = WWFR.map_mc(model, sample, reduce!, 20)

println("WhereTheWaterFlowsRandomly simple example (WWF + uncertain source)")
println("samples: ", aggr.n_samples[])
println("mean catchment discharge: ", mean(aggr.catchment_fluxes[1]))
println("standard deviation of catchment discharge: ", std(aggr.catchment_fluxes[1]))
println("min/max catchment discharge: ", Float64.(extrema(aggr.catchment_fluxes[1])))
