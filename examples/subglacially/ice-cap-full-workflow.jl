# examples/subglacially/ice-cap-full-workflow.jl
#
# Full subglacial routing workflow for an ice cap with multiple outlet glaciers.
#
# Demonstrates on a synthetic Vatnajökull-like geometry:
#   1. Deterministic subglacial routing with per-outlet catchment diagnostics.
#   2. The two-step ctch_sinks pattern: discover actual routing sinks first,
#      then partition them into outlet groups.
#   3. Monte Carlo bed-uncertainty propagation and stochastic catchment-flux
#      analysis.
#
# Run from the examples/ environment:
#   include("subglacially/ice-cap-full-workflow.jl")

using WhereTheWaterFlows, CairoMakie, Statistics, Random

const WWFS = WhereTheWaterFlows.Subglacially
const WWFR = WhereTheWaterFlows.Randomly

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Synthetic Vatnajökull-like data
#
# Vatnajökull is roughly 150 km × 100 km.  We use a 100×80 grid at 1500 m
# spacing.  All coordinates are in metres.
# ─────────────────────────────────────────────────────────────────────────────
Random.seed!(43)

nx, ny = 100, 80
dx     = 1_500.0              # grid spacing [m]
x      = collect(range(-nx÷2*dx, step=dx, length=nx))
y      = collect(range(-ny÷2*dx, step=dx, length=ny))

# Dome-shaped ice surface, asymmetric (peak offset from centre, like Vatnajökull)
cx = cy = 0.0
ccx = 3.5e4   # trough x-offset [m]
ccy = 3.5e4   # trough y-offset [m]

R2 = 55_000.0^2               # dome radius squared [m²]

surfdem = Float32[
    1800.0 * max(1.0 - ((xi - cx)^2 + (yi - cy)^2) / R2, 0.0)^0.4
    for xi in x, yi in y
]

# Bed topography: lower than surface, with two subglacial troughs leading
# toward the south and east margins.
beddem = Float32[
    ( 400.0 * max(1.0 - ((xi - cx)^2 + (yi - cy)^2) / 72_000.0^2, 0.0)^0.3
      -  80.0 * exp(-((xi - 0.55*ccx)^2 + (yi - 1.30*ccy)^2) / 10_000.0^2)
      -  80.0 * exp(-((xi - 0.55*ccx)^2 + (yi + 1.30*ccy)^2) / 10_000.0^2)
      -  60.0 * exp(-((xi + 1.40*ccx)^2 + (yi - 0.55*ccy)^2) / 12_000.0^2) )
    for xi in x, yi in y ] .+ 10.0 .* Float32.(randn(nx, ny))

# Enforce surface ≥ bed everywhere
surfdem .= max.(surfdem, beddem .+ 1.0f0)

# Glacier mask: cells with ice thickness > 50 m
thick   = surfdem .- beddem
glacier = thick .> 50.0f0

# Replace off-glacier bed with surface elevation so that the Shreve hydraulic
# potential φ = f·H·(ρᵢ/ρ_w) + z_b stays finite everywhere.  The routing
# mask then restricts actual routing to glacier cells only.
beddem_clean = copy(beddem)
beddem_clean[.!glacier] .= surfdem[.!glacier]

# Basal melt source: uniform 5 mm/day converted to m/s.
# `source` is per unit area [m/s]; the routing accumulates source × dx² [m³/s].
source = fill(Float32(5e-3 / 86_400), nx, ny)

println("Synthetic Vatnajökull grid: $(nx)×$(ny), dx=$(dx) m")
println("Ice-covered cells: $(sum(glacier)) / $(nx*ny) " *
        "($(round(100*mean(glacier), digits=1)) %)")

# ── Plot 1: Input fields ──────────────────────────────────────────────────────
surfdem_plt = copy(surfdem); surfdem_plt[.!glacier] .= NaN
beddem_plt  = copy(beddem);  beddem_plt[.!glacier]  .= NaN
thick_plt   = copy(thick);   thick_plt[.!glacier]   .= NaN

fig1 = Figure(size=(1100, 320))
heatmap(fig1[1, 1], x, y, surfdem_plt;
        colormap=:ice,     axis=(title="Surface elevation [m]",
                                 xlabel="x [m]", ylabel="y [m]"))
heatmap(fig1[1, 2], x, y, beddem_plt;
        colormap=:terrain, axis=(title="Bed elevation [m]",
                                 xlabel="x [m]"))
heatmap(fig1[1, 3], x, y, thick_plt;
        colormap=:viridis, axis=(title="Ice thickness [m]",
                                 xlabel="x [m]"))
fig1
# save("ice-cap_1_input.png", fig1)

# ─────────────────────────────────────────────────────────────────────────────
# 2.  First routing pass — discover actual routing sinks
#
# `ctch_sinks` must be built from the *actual* routing sinks (cells where flow
# physically exits the active domain), not from an arbitrary set of boundary
# indices.  Off-glacier cells are barriers, not sinks: pre-building groups from
# raw CartesianIndices would silently include zero-flux cells.
#
# Two-step approach:
#   a) Run routing once without ctch_sinks to obtain out.routing.sinks.
#   b) Partition those sinks by position into outlet groups.
#   c) Re-run with the partitioned ctch_sinks (or pass directly to map_mc).
#
# For real datasets, WWFS.catchment_sinks automates step (b) from polygon
# outlines of each outlet glacier.
# ─────────────────────────────────────────────────────────────────────────────
println("\n── First pass: discovering active sinks ─────────────────────────────")

out_first = WWFS.waterflows_subglacial(
    surfdem, beddem_clean, dx,
    1.0f0, source, glacier;
    gamma       = WWFS.GAMMA,
    drain_pits  = true,
    bnd_as_sink = true,
    nan_as_sink = true,
)

all_sinks = out_first.routing.sinks
println("Total active sinks discovered: ", length(all_sinks))

# Partition sinks into four quadrant-based outlet sectors.
# ci[1] is the x-index, ci[2] is the y-index.
sw_outlet = filter(ci -> ci[1] <= nx÷2 && ci[2] <= ny÷2, all_sinks)
se_outlet = filter(ci -> ci[1] >  nx÷2 && ci[2] <= ny÷2, all_sinks)
nw_outlet = filter(ci -> ci[1] <= nx÷2 && ci[2] >  ny÷2, all_sinks)
ne_outlet = filter(ci -> ci[1] >  nx÷2 && ci[2] >  ny÷2, all_sinks)

ctch_sinks   = [sw_outlet, se_outlet, nw_outlet, ne_outlet]
outlet_names = ["SW (Skeiðará-like)", "SE", "NW (terminus)", "NE (Jökulsá-like)"]

println("Active sinks per outlet:")
for (name, grp) in zip(outlet_names, ctch_sinks)
    println("  $(rpad(name, 24)) $(length(grp))")
end

# Verify the partition is exhaustive (non-overlapping groups covering all sinks)
@assert length(all_sinks) == sum(length.(ctch_sinks))

# ─────────────────────────────────────────────────────────────────────────────
# 3.  Deterministic subglacial routing (with outlet groups)
# ─────────────────────────────────────────────────────────────────────────────
println("\n── Deterministic subglacial routing ─────────────────────────────────")

out = WWFS.waterflows_subglacial(
    surfdem, beddem_clean, dx,
    1.0f0,      # floatfrac: full flotation everywhere
    source,     # basal melt [m/s]
    glacier;    # routing mask
    gamma       = WWFS.GAMMA,   # Röthlisberger constant ≈ −0.31
    ctch_sinks  = ctch_sinks,
    drain_pits  = true,
    bnd_as_sink = true,
    nan_as_sink = true,
)

(; sinks, pits) = out.routing

println("Total outlet sinks:         ", length(sinks))
println("Unresolved interior pits:   ", length(pits))
println("Supercooled cells:          ", sum(out.pressmelt.sc_locs))
println("Max discharge [m³/s]:       ",
        round(maximum(out.routing.area.total[glacier]), sigdigits=3))
println("Max lake depth (fixed) [m]: ",
        round(maximum(out.lakes.depth_fixed_surface[glacier]), sigdigits=3))

# out.sink_catchments.fluxes is a named tuple (total, dissipation, pressmelt),
# each a Vector with one scalar per ctch_sinks group.
println("\nPer-outlet flux (deterministic) [m³/s]:")
for (name, q) in zip(outlet_names, out.sink_catchments.fluxes.total)
    println("  $(rpad(name, 24)) $(round(q, sigdigits=3))")
end

# ── Plot 2a: Discharge (log₁₀) with sink locations ───────────────────────────
fig2a = plt_area(x, y, out.routing.area.total; sinks=out.routing.sinks)
fig2a
# save("ice-cap_2a_discharge.png", fig2a)

# ── Plot 2b: Catchment map ────────────────────────────────────────────────────
fig2b = plt_catchments(x, y, out.routing.c; minsize=1)
# save("ice-cap_2b_catchments.png", fig2b)

# ── Plot 3: Subglacial diagnostics ───────────────────────────────────────────
# Hydraulic potential, lake depth, and supercooling mask.
phi_plt  = copy(out.routing.phi);               phi_plt[.!glacier]  .= NaN
lake_plt = copy(out.lakes.depth_fixed_surface); lake_plt[.!glacier] .= NaN
sc_plt   = Float32.(out.pressmelt.sc_locs);     sc_plt[.!glacier]   .= NaN

fig3 = Figure(size=(1100, 320))
heatmap(fig3[1, 1], x, y, phi_plt;
        colormap=:curl,  axis=(title="Shreve hydraulic potential φ [m]",
                               xlabel="x [m]", ylabel="y [m]"))
heatmap(fig3[1, 2], x, y, lake_plt;
        colormap=:blues, axis=(title="Lake depth — fixed surface [m]",
                               xlabel="x [m]"))
heatmap(fig3[1, 3], x, y, sc_plt;
        colormap=:hot,   axis=(title="Supercooled cells",
                               xlabel="x [m]"))
fig3
# save("ice-cap_3_subglacial.png", fig3)

# ── Plot 3b: Lake depth via fill_dem ────────────────────────────────────
phi_filled = fill_dem(out.routing.phi, out.routing.sinks, out.routing.dir)
fig3b = heatmap(x, y, phi_filled .- out.routing.phi)
# save("ice-cap_3b_lakedepth.png", fig3b)

# ─────────────────────────────────────────────────────────────────────────────
# 4.  Monte Carlo bed-uncertainty propagation
#
# Propagate 10 % relative bed uncertainty with a 3 km spatial correlation
# length (roughly one ice thickness for this synthetic glacier) to obtain a
# stochastic distribution of per-outlet catchment fluxes.
#
# Workflow:
#   a) Describe each uncertain field with WWFR.Uncertainty.
#   b) Build model / sample / reduce! with WWFR.make_fns_subglacial.
#   c) Run the Monte Carlo loop with WWFR.map_mc.
# ─────────────────────────────────────────────────────────────────────────────
println("\n── Monte Carlo bed-uncertainty propagation ───────────────────────────")

surfdem_uc   = WWFR.Uncertainty()   # surface treated as known exactly
beddem_uc    = WWFR.Uncertainty(reluc=0.10, correlation_length=3_000.0)
floatfrac_uc = WWFR.Uncertainty()   # flotation fraction treated as exact
source_uc    = WWFR.Uncertainty()   # melt source treated as exact

# make_fns_subglacial requires floatfrac as an array (not a scalar).
floatfrac_arr = fill(Float32(1.0), nx, ny)

model, sample, reduce! = WWFR.make_fns_subglacial(
    dx,
    surfdem,       surfdem_uc,
    beddem_clean,  beddem_uc,
    floatfrac_arr, floatfrac_uc,
    source,        source_uc,
    ctch_sinks;
    mask  = glacier,
    gamma = WWFS.GAMMA,
)

# 20 samples is deliberately small for quick execution.
# Use 50–200 samples for stable statistics in production.
n_samples = 20
Random.seed!(123)
aggr = WWFR.map_mc(model, sample, reduce!, n_samples)

println("Samples completed:            ", aggr.n_samples[])
println("Mean max discharge [m³/s]:    ",
        round(maximum(aggr.areas_total[glacier]), sigdigits=3))
println("Max supercooling frequency:   ",
        round(maximum(aggr.sc_locs[glacier]),     sigdigits=3))
println("Max lake-occurrence fraction: ",
        round(maximum(aggr.lake_mask_fixed_surface[glacier]), sigdigits=3))

# ── Plot 4: MC mean spatial fields ───────────────────────────────────────────
mc_area_plt = copy(aggr.areas_total);             mc_area_plt[.!glacier] .= NaN
mc_sc_plt   = copy(aggr.sc_locs);                mc_sc_plt[.!glacier]   .= NaN
mc_lake_plt = copy(aggr.lake_mask_fixed_surface); mc_lake_plt[.!glacier] .= NaN

fig4 = Figure(size=(1100, 320))
heatmap(fig4[1, 1], x, y, log10.(max.(mc_area_plt, 1f-10));
        colormap=:viridis, axis=(title="MC mean discharge log₁₀ [m³/s]",
                                  xlabel="x [m]", ylabel="y [m]"))
heatmap(fig4[1, 2], x, y, mc_sc_plt;
        colormap=:inferno, colorrange=(0, 1),
        axis=(title="Supercooling frequency (0–1)", xlabel="x [m]"))
heatmap(fig4[1, 3], x, y, mc_lake_plt;
        colormap=:Blues,   colorrange=(0, 1),
        axis=(title="Lake-occurrence fraction (0–1)", xlabel="x [m]"))
fig4
# save("ice-cap_4_mc_fields.png", fig4)

# Stochastic catchment-flux analysis.
#
# aggr.catchment_fluxes for subglacial MC is a named tuple
# (total, dissipation, pressmelt), mirroring sink_catchments.fluxes from the
# deterministic run.  aggr.catchment_fluxes.total[i] is a Vector{Float32} of
# length n_samples containing the total flux at outlet i in each realization.
# (The subaerial counterpart uses aggr.catchment_fluxes[i] directly, with no
# .total field.)
println("\nPer-outlet catchment flux statistics [m³/s] ($(n_samples) samples):")
println("  $(rpad("Outlet", 26)) mean         std          CV")
for (i, name) in enumerate(outlet_names)
    fluxes = aggr.catchment_fluxes.total[i]
    μ = mean(fluxes)
    σ = std(fluxes)
    println("  $(rpad(name, 26)) " *
            "$(rpad(round.(μ, sigdigits=3), 13))" *
            "$(rpad(round.(σ, sigdigits=3), 13))" *
            "$(round.(σ/μ, sigdigits=2))")
end

# ── Plot 5: Per-outlet catchment flux distributions ──────────────────────────
# Bar chart (mean ± 1 σ) with individual sample points overlaid.
fig5 = Figure(size=(700, 400))
ax5  = Axis(fig5[1, 1];
            title    = "Per-outlet discharge — MC spread ($(n_samples) samples)",
            xlabel   = "Outlet",
            ylabel   = "Total flux [m³/s]",
            xticks   = (1:length(outlet_names), outlet_names),
            xticklabelrotation = π/6)

colors = Makie.wong_colors()
for (i, name) in enumerate(outlet_names)
    fluxes = aggr.catchment_fluxes.total[i]
    μ = mean(fluxes)
    σ = std(fluxes)
    barplot!(ax5, [i], [μ]; color=(colors[i], 0.6), width=0.6)
    errorbars!(ax5, [i], [μ], [σ]; whiskerwidth=10, color=:black)
    scatter!(ax5, fill(i, length(fluxes)), fluxes;
             color=(colors[i], 0.9), markersize=6)
end
fig5
# save("ice-cap_5_outlet_fluxes.png", fig5)

# ── Plot 6: Per-outlet flux histograms (2×2) ─────────────────────────────────
fig6 = Figure(size=(700, 600))
for (i, name) in enumerate(outlet_names)
    row, col = fldmod1(i, 2)
    ax = Axis(fig6[row, col];
              title  = name,
              xlabel = "Total flux [m³/s]",
              ylabel = "Count")
    hist!(ax, aggr.catchment_fluxes.total[i]; color=(Makie.wong_colors()[i], 0.7))
end
fig6
# save("ice-cap_6_outlet_histograms.png", fig6)
