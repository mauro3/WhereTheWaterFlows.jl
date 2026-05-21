"""
This script calculates the channel deflection angle which maximizes melt.

It does this for the continuous case as well as the "quantized" D8 routing case.

Variables:
- R: ratio between |-∇ϕ_m| and |-∇ϕ|, i.e. between bed slope and hydraulic slope
- lambda: angle between -∇ϕ and -∇ϕ_m
- kappa: angle between -∇ϕ and the flow direction
"""

using WhereTheWaterFlows
const WWFS = WhereTheWaterFlows.Subglacially
const WWF = WhereTheWaterFlows
using CairoMakie
import Roots, Optim

"Röthlisberger constant"
const gamma = WWFS.GAMMA
"Supercooling threshold"
sc_threshold = (1+gamma) / gamma

"""
Function to be maximised to get angle kappa with max melt.
"""
function L(kappa, R, lambda)
    beta = [3/2, 2][1]
    if kappa<0 || kappa>pi/2
        return NaN
    end
    (cos(kappa)^(beta-1)) * ( (1+gamma)*cos(kappa) - gamma * cos(lambda-kappa) * R)
end

"""
Return kappa (angle between -∇ϕ and Q) which maximizes melt.
"""
function kappa(R, lambda)
    if lambda<0
        fac = -1
        lambda = -lambda
    else
        fac = 1
    end
    sol = Optim.optimize(kappa -> -L(kappa[1], R, lambda), [1e-3], Optim.LBFGS())
    if Optim.converged(sol)
        return fac * Optim.minimizer(sol)[1]
    else
        return NaN
    end
end

"""Ratio of melt in a "wonkey" channel vs one on a flat bed (i.e. co-linear with -∇ϕ)"""
meltratio(kappa_max, R, lambda) = L(kappa_max, R, lambda) / (1+gamma)
meltratio(R, lambda) = meltratio(kappa(R, lambda), R, lambda)

"""
Finds locus where the meltratio==mu as a function of R.
"""
function locus(mu, R)
    bracket = (0.0,pi)
    f = lambda -> meltratio(R, lambda) - mu
    if sign(f(bracket[1])) == sign(f(bracket[2]))
        return NaN
    else
        return Roots.find_zero(f, bracket)
    end
end


# Plotting etc
Rs = 0:0.02:8.5
lambdas = 0:0.02:pi
kappas = kappa.(Rs, lambdas');
mr = meltratio.(kappas, Rs, lambdas');

mr_c = copy(mr)
mr_c[mr.<0.05] .= NaN
lambdas1 = locus.(1.0, Rs[1:end])
lambdas01 = locus.(0.1, Rs[1:end])

# Plot deflection angle kappa and melt ratio m in terms of R and lambda
fig = Figure();
ax1 = Axis(fig[1,1], xlabel = "λ (°)", ylabel = "R", limits = (0, 180, 0, 8.5))
cf = contourf!(ax1, lambdas.*180/pi, Rs, kappas'.*180/pi,levels = 0:15/4:90, colormap=:blues)
Colorbar(fig[1,2], cf,ticks = 0:15:90,label = "κ (°)")
lines!(lambdas1.*180/pi, Rs, label="m=1", color=:black, linestyle=:solid)
lines!(lambdas01.*180/pi, Rs, label="m=0.1", color=:black, linestyle=:dash)
axislegend(position=:lt)
ax2 = Axis(fig[1,3], xlabel = "λ (°)", ylabel = "R", limits = (0, 180, 0, 8.5))
# https://discourse.julialang.org/t/custom-divergent-colormap-with-unequal-sides/115484/7
cm = Makie.diverging_palette(0, 230; c=0.6, mid=1/5)
p2 = contourf!(ax2, lambdas.*180/pi, Rs, mr', levels=0:0.25:5, colormap=cm) #, xlabel="λ (°)", ylabel="R",
Colorbar(fig[1,4], p2, ticks = 0:0.5:5,label = "m ()")
lines!(lambdas1.*180/pi, Rs, label="m=1", color=:black, linestyle=:solid)
lines!(lambdas01.*180/pi, Rs, label="m=0.1", color=:black, linestyle=:dash)
axislegend(position=:lt)
Label(fig[1,1,TopLeft()], "A", font=:bold)
Label(fig[1,3,TopLeft()], "B", font=:bold)
save(joinpath(@__DIR__, "deflection.png"), fig)

## Look at lookup table for the D8 case
fig = Figure();
ax1 = Axis(fig[1,1], xlabel = "R ()", ylabel = "κ (°)", limits = (0, 8.5, 0, 95))
ratio =0:0.01:8.5
cs = [:red, :magenta, :green]
for (i, lambda) = enumerate((3:-1:1).*pi/4)
    kappa_ = kappa.(ratio, lambda)
    kappa__ = round.(Int, kappa_ / (pi/4))*pi/4
    kappa__ = kappa__ .+ pi/2/180 * (1-i) * 1.0
    lines!(ratio, rad2deg.(kappa__), color=cs[i], linestyle=:solid,
            label="λ=$(rad2deg(lambda)), D8-quantized")
    lines!(ratio, rad2deg.(kappa_), color=cs[i], linestyle=:dash,
            label="                   continuous")
end
axislegend(labelsize=9, position=:lt, rowgap=0.1, height=130)
# cg = cgrad(:Blues_9)
ax2 = Axis(fig[1,2], xlabel = "λ (°)", ylabel = "R", limits = (0, 180, 0, 8.5))
cf = contourf!(ax2, lambdas.*180/pi, Rs, kappas'.*180/pi,levels = 0:15/4:90, colormap=:blues)
# cg = cgrad(:bone)[1:10]
# contours corresponding to where deflection angles change to 45 and 90deg, respectively
contour!(ax2, lambdas.*180/pi, Rs, kappas'.*180/pi, levels=[45/2, 45*3/2], linestyle=:dash, color=:black)
Colorbar(fig[1,3], cf,ticks = 0:15:90,label = "κ (°)")
# vertical lines corresponding to D8 directions
lines!(ax2, Rs.*0 .+ 45, Rs, label="", linestyle=:dash, color=cs[3])
lines!(ax2, Rs.*0 .+ 90, Rs, label="", linestyle=:dash, color=cs[2])
lines!(ax2, Rs.*0 .+ 135, Rs, label="", linestyle=:dash, color=cs[1])
Label(fig[1,1,TopLeft()], "A", font=:bold)
Label(fig[1,2,TopLeft()], "B", font=:bold)
save(joinpath(@__DIR__, "deflection-D8.png"), fig)
