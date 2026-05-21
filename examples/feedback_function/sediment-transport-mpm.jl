using CairoMakie
using Random; Random.seed!(42)
using WhereTheWaterFlows

const WWF = WhereTheWaterFlows

# Small synthetic catchment with a clear downslope direction.
n = 100
dx = 50.0
x = y = collect(0:dx:(n - 1) * dx)
dem = 300 .- 0.02.*x .+ 0.003.*y' .+ 0.8.*sin.(x./500) .+ 0.3.*randn(n, n)

# Inputs routed by WWF:
# - water discharge source [m^3/s per cell]
# - sediment source [arbitrary extensive unit per cell]
water_source = fill(0.02, n, n)
sediment_source = fill(5e-4, n, n)
cellarea = (water_source, sediment_source)

# Simple MPM-like closure with width scaling W ~ Q^b.
const W0 = 1.0
const b_width = 0.5
const Qref = 1.0
const theta_c = 0.1
const k_theta = 45.0
const alpha = 2.5

function mpm_capacity(Q, slope)
    W = max(W0 * (Q / Qref)^b_width, 0.5)
    theta = k_theta * slope
    qb = theta > theta_c ? alpha * (theta - theta_c)^(3/2) : 0.0
    return qb * W
end

function sediment_feedback(uparea, ij, dir)
    Q, Qs_in = uparea
    d = dir[ij]

    d >= WWF.SINK && return (Q, Qs_in)

    ds = iseven(d) ? dx : sqrt(2.0) * dx
    ij2 = ij + WWF.dir2ind(d)
    slope = max((dem[ij] - dem[ij2]) / ds, 0.0)

    Qs_cap = mpm_capacity(max(Q, 1e-8), slope)
    Qs_out = min(Qs_in, Qs_cap)
    return (Q, Qs_out)
end

out_no_feedback = WWF.waterflows(dem, cellarea; drain_pits=true)
out_mpm = WWF.waterflows(dem, cellarea; drain_pits=true, feedback_fn=sediment_feedback)

Q = out_mpm.area[1]
Qs = out_mpm.area[2]

sed_in = sum(sediment_source)
sed_out_no_feedback = sum(out_no_feedback.area[2][out_no_feedback.sinks])
sed_out_mpm = sum(Qs[out_mpm.sinks])

println("Simple MPM-like sediment transport")
println("total sediment input: ", sed_in)
println("sediment outflow without capacity limit: ", sed_out_no_feedback)
println("sediment outflow with MPM-like capacity: ", sed_out_mpm)

fig = Figure(size=(980, 360))

ax1 = Axis(fig[1, 1], title="Water discharge Q")
hm1 = heatmap!(ax1, x, x, log10.(max.(Q, 1e-8)))
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[1, 3], title="Sediment Qs (no feedback)")
hm2 = heatmap!(ax2, x, x, log10.(max.(out_no_feedback.area[2], 1e-12)))
Colorbar(fig[1, 4], hm2)

ax3 = Axis(fig[1, 5], title="Sediment Qs (MPM-like)")
hm3 = heatmap!(ax3, x, x, log10.(max.(Qs, 1e-12)))
Colorbar(fig[1, 6], hm3)

display(fig)
