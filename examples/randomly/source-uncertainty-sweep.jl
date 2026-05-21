using WhereTheWaterFlows
using Random; Random.seed!(42)
using Statistics
const WWFR = WhereTheWaterFlows.Randomly

function build_dem(n)
    t = range(-pi, pi, length=n)
    dem = sin.(t) .* cos.(t')
    return dem
end

function run_case(absuc, reluc, ondem; n_samples=24, n=100, dx=100.0)
    dem = build_dem(n)
    source = fill(1.0 / dx^2, size(dem))
    dem_uc = WWFR.Uncertainty(absuc=absuc*ondem, reluc=reluc*ondem, correlation_length=10dx)
    source_uc = WWFR.Uncertainty(absuc=ondem * absuc / dx^2,
                                 reluc=ondem * reluc,
                                 correlation_length=10dx)

    ctch_sinks = [CartesianIndices((2:2, 2:n-1))[:]]

    model, sample, reduce! = WWFR.make_fns_subaerial(dx,
                                                     dem, dem_uc,
                                                     source, source_uc,
                                                     ctch_sinks)

    aggr = WWFR.map_mc(model, sample, reduce!, n_samples, progressmeter=false)

    max_area = maximum(aggr.areas_total)
    max_slen = maximum(aggr.stream_length)
    mean_sink_flux = mean(aggr.catchment_fluxes[1])
    std_sink_flux = std(aggr.catchment_fluxes[1])
    # display(plt_area(1:n,1:n,aggr.areas_total))

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
    out = run_case(absuc, reluc, false)
    println(lpad(round(out.absuc, digits=2), 6),
            lpad(round(out.reluc, digits=2), 6),
            lpad(round(out.max_area, digits=2), 11),
            lpad(round(out.max_slen, digits=2), 11),
            lpad(round(out.mean_sink_flux, digits=0), 17),
            lpad(round(out.std_sink_flux, digits=0), 16))
end

println("WWFR source-uncertainty sweep (WWF routing)")
println("  uncertain DEM, fixed source term")
println()
println(" absuc reluc   max_area   max_slen   mean_sink_flux   std_sink_flux")

for (absuc, reluc) in cases
    out = run_case(absuc, reluc, true)
    println(lpad(round(out.absuc, digits=2), 6),
            lpad(round(out.reluc, digits=2), 6),
            lpad(round(out.max_area, digits=2), 11),
            lpad(round(out.max_slen, digits=2), 11),
            lpad(round(out.mean_sink_flux, digits=0), 17),
            lpad(round(out.std_sink_flux, digits=0), 16))
end
