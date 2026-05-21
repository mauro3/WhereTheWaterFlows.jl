# Examples

Run examples from this folder:

```julia
julia --project
include("wwf-simple.jl")
```

## Quick Start (top-level)

- `wwf-simple.jl`: core deterministic routing with `WhereTheWaterFlows`
- `wwfs-simple.jl`: subglacial deterministic routing with `WhereTheWaterFlows.Subglacially`
- `wwfr-simple.jl`: stochastic routing with `WhereTheWaterFlows.Randomly` (uncertain source, fixed DEM)

## Topical Folders

- `core/`: additional WWF workflows and plotting demos
- `subglacially/`: richer subglacial scenarios
- `randomly/`: advanced stochastic workflows
- `theory/`: theory and diagnostics
- `data/`: example data and provenance scripts

## Suggested order

1. `wwf-simple.jl`
2. `wwfs-simple.jl`
3. `wwfr-simple.jl`
4. Then browse a topical folder of interest.
