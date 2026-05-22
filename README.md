# WhereTheWaterFlows

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://mauro3.github.io/WhereTheWaterFlows.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://mauro3.github.io/WhereTheWaterFlows.jl/dev)

[![Build Status](https://github.com/mauro3/WhereTheWaterFlows.jl/workflows/CI/badge.svg)](https://github.com/mauro3/WhereTheWaterFlows.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mauro3/WhereTheWaterFlows.jl?svg=true)](https://ci.appveyor.com/project/mauro3/WhereTheWaterFlows-jl)
[![Coverage](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl)
[![DOI](https://zenodo.org/badge/218504028.svg)](https://doi.org/10.5281/zenodo.7086860)

This package allows to calculate water flow paths on digital elevation models (DEMs).

This package implements the D8 flow routing algorithm [1] as well as a
basin-filling algorithm, also by [1].  In its implementation it uses a
O(n), recursive algorithm similar as in [2].  Due to recursion it can run
into a stackoverflow error on very large DEMs.

This code is reasonably fast: flow routing on a DEM of Antarctica of
about 2e8 points (14000x14000) with 150000 depressions takes about 30s
on a Ryzen 4750U laptop.

WhereTheWaterFlows has been used in glaciological contexts [3,4,5,6] but could be useful in other settings as well.

# Example

Example of an upslope area calculation on a synthetic digital elevation model:

```julia
using WhereTheWaterFlows, CairoMakie

# build a small synthetic DEM
n = 200
xs = range(-π, π, length=n)
dem = sin.(xs) .* cos.(xs') .+ 0.05 * rand(n,n)

out = waterflows(dem)
plt_area(xs, xs, out.area)
```
![Upslope area](https://user-images.githubusercontent.com/4098145/67853636-e319b880-fb06-11e9-933d-9f55ace99ce1.png)

For details, see the [documentation](https://mauro3.github.io/WhereTheWaterFlows.jl).

# References
[1] O’Callaghan, J. and Mark, D.: The extraction of drainage networks
    from digital elevation data, Comput. Vision Graph., 28, 323–344,
    1984. [download via google scholar](https://scholar.google.ch/scholar?hl=en&as_sdt=0%2C5&q=The+extraction+of+drainage+networks+from+digital+elevation+data&btnG=)

[2] Braun, J. and Willett, S.D.: A very efficient O(n), implicit and
    parallel method to solve the stream power equation governing
    fluvial incision and landscape evolution, Geomorphology, 2013.
    [link](https://doi.org/10.1016/j.geomorph.2012.10.008)

[3] Malczyk et al., Constraints on Subglacial Melt Fluxes from Observations of Active Subglacial Lake Recharge, Journal of Glaciology, 2023. [link](https://doi.org/10.1017/jog.2023.70)

[4] Delaney et al., Modeling the Spatially Distributed Nature of Subglacial Sediment Transport and Erosion, Earth Surface Dynamics, 2023. [link](https://doi.org/10.5194/esurf-11-663-2023)

[5] Ogier et al., Definition, formation and rupture mechanisms of water pockets in alpine glaciers: Insights from an updated inventory for the Swiss Alps, Journal of Glaciology, 2025. [link](https://doi.org/10.1017/jog.2025.43)

[6] Ogier et al., Potential Glacier Contributions to the 2024 La Bérarde Flood, EGUsphere pre-print, 2026. [link](https://doi.org/10.5194/egusphere-2026-466)
