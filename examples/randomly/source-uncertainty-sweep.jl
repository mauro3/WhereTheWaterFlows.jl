using WhereTheWaterFlows
using Random
using Statistics

const WWFR = WhereTheWaterFlows.Randomly

Random.seed!(42)

function build_dem(n)
    xs = range(-pi, pi, length=n)
    dem = sin.(xs) .* cos.(xs')
    mask = trues(size(dem))
    return xs, dem, mask
end

function run_case(absuc, reluc; n_samples=24, n=100, dx=100.0)
    _, dem, mask = build_dem(n)
    source = fill(1.0 / dx^2, size(dem))
    dem_uc = WWFR.Uncertainty(absuc=0.0, reluc=0.0, correlation_length=dx)
    source_uc = WWFR.Uncertainty(absuc=absuc / dx^2,
                                 reluc=reluc,
                                 correlation_length=8.0 * dx)

    ctch_sinks = [CartesianIndices((2:2, 2:n-1))[:]]

    model, sample, reduce! = WWFR.make_fns_subaerial(dx,
                                                     dem, dem_uc,
                                                     source, source_uc,
                                                     ctch_sinks,
                                                     mask)

    aggr = WWFR.map_mc(model, sample, reduce!, n_samples)

    max_area = maximum(aggr.areas_total)
    max_slen = maximum(aggr.stream_length)
    mean_sink_flux = mean(aggr.catchment_fluxes[1])
    std_sink_flux = std(aggr.catchment_fluxes[1])

    return (; absuc, reluc, max_area, max_slen, mean_sink_flux, std_sink_flux)
end

cases = [
    (0.00, 0.00),
    (0.10, 0.00),
    (0.20, 0.00),
    (0.00, 0.10),
    (0.00, 0.20),
    (0.20, 0.20),
]

println("WWFR source-uncertainty sweep (WWF routing)")
println("  fixed DEM, uncertainty only in source term")
println()
println(" absuc reluc   max_area   max_slen   mean_sink_flux   std_sink_flux")

for (absuc, reluc) in cases
    out = run_case(absuc, reluc)
    println(lpad(round(out.absuc, digits=2), 6),
            lpad(round(out.reluc, digits=2), 6),
            lpad(round(out.max_area, digits=2), 11),
            lpad(round(out.max_slen, digits=2), 11),
            lpad(round(out.mean_sink_flux, digits=4), 17),
            lpad(round(out.std_sink_flux, digits=4), 16))
end
