# Examples

Run examples from this folder with:

```julia
julia --project
using Pkg; Pkg.up() # only run the very first time to install packages
include("wwf-simple.jl") # or any other script
```

Introduction scripts (no plotting):

- `wwf-simple.jl`: core deterministic routing with `WhereTheWaterFlows`
- `wwfs-simple.jl`: subglacial deterministic routing with `WhereTheWaterFlows.Subglacially`
- `wwfr-simple.jl`: stochastic routing with `WhereTheWaterFlows.Randomly` (uncertain source, fixed DEM)

Folders with more in-depth scripts:

- `core/`: additional WWF workflows and plotting demos
- `subglacially/`: richer subglacial scenarios
- `randomly/`: advanced stochastic workflows
- `theory/`: theory and diagnostics
- `data/`: example data and provenance scripts

