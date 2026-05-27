# WhereTheWaterFlows

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://mauro3.github.io/WhereTheWaterFlows.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://mauro3.github.io/WhereTheWaterFlows.jl/dev)

[![Build Status](https://github.com/mauro3/WhereTheWaterFlows.jl/workflows/CI/badge.svg)](https://github.com/mauro3/WhereTheWaterFlows.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mauro3/WhereTheWaterFlows.jl?svg=true)](https://ci.appveyor.com/project/mauro3/WhereTheWaterFlows-jl)
[![Coverage](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mauro3/WhereTheWaterFlows.jl)
[![DOI](https://zenodo.org/badge/218504028.svg)](https://doi.org/10.5281/zenodo.7086860)

This package calculates water flow paths on digital elevation models (DEMs).

It implements the D8 flow routing algorithm combined with a
breach-type basin-filling algorithm as descibed by [1].  Its implementation uses a
O(n), recursive algorithm [2].

This model is reasonably fast: flow routing on a DEM of Antarctica of
about 2e8 points (14000x14000) with 150000 depressions takes about 30s
on a Ryzen 4750U laptop. Due to recursion it can run
into a stackoverflow error on very large DEMs.

WhereTheWaterFlows has been used in glaciological contexts [3,4,5,6] but could be useful in other settings as well.

# Example

Example of an upslope area calculation on a synthetic digital elevation model:

```julia
using WhereTheWaterFlows, CairoMakie

# build a small synthetic DEM
n = 200
x = y = range(-π, π, length=n)
dem = sin.(x).*cos.(y') .+ 0.05*rand(n,n)

out = waterflows(dem)
plt_area(x, y, out.area)
```
![Upslope area](https://private-user-images.githubusercontent.com/4098145/598786029-5bfd3617-07c0-4d4c-b5d8-ca9e91a4d956.png?jwt=eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTUiLCJleHAiOjE3Nzk4ODQyMzQsIm5iZiI6MTc3OTg4MzkzNCwicGF0aCI6Ii80MDk4MTQ1LzU5ODc4NjAyOS01YmZkMzYxNy0wN2MwLTRkNGMtYjVkOC1jYTllOTFhNGQ5NTYucG5nP1gtQW16LUFsZ29yaXRobT1BV1M0LUhNQUMtU0hBMjU2JlgtQW16LUNyZWRlbnRpYWw9QUtJQVZDT0RZTFNBNTNQUUs0WkElMkYyMDI2MDUyNyUyRnVzLWVhc3QtMSUyRnMzJTJGYXdzNF9yZXF1ZXN0JlgtQW16LURhdGU9MjAyNjA1MjdUMTIxMjE0WiZYLUFtei1FeHBpcmVzPTMwMCZYLUFtei1TaWduYXR1cmU9YjU5MDc2NmUwZjMzNjE4ZDhhZmUxZWVlZGZhNTEwYjNkYjFhMWI3NDhlYWViM2U1YTQ4YzZlOGY0NTA3NWM4YyZYLUFtei1TaWduZWRIZWFkZXJzPWhvc3QmcmVzcG9uc2UtY29udGVudC10eXBlPWltYWdlJTJGcG5nIn0.LRG_LqUv7KJ7o0KXNY13OZGg4KPX1UntBH3wDYREbH0)

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
