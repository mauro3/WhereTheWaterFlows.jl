# Examples

Run examples from this folder with:

```julia
julia --project
using Pkg; Pkg.up() # only run the very first time to install packages
include("wwf-simple.jl") # or any other script
```

Run all examples in the `examples` environment with:

```julia
include("run-all-examples.jl")
```

Introduction scripts (no plotting):

- `wwf-simple.jl`: core deterministic routing with `WhereTheWaterFlows`
- `wwfs-simple.jl`: subglacial deterministic routing with `WhereTheWaterFlows.Subglacially`
- `wwfr-simple.jl`: stochastic routing with `WhereTheWaterFlows.Randomly` (uncertain source, fixed DEM)

Folders with more in-depth scripts:

- `core/`: additional WWF workflows and plotting demos
- `feedback_function/`: shows how the `feedback_fn` of `waterflows` can be used for grid-based computations
  - includes melt-production and MPM-like sediment transport examples
- `subglacially/`: subglacial scenarios using WWFS
  - includes `subglacially/ice-cap-full-workflow.jl`: end-to-end deterministic
    plus Monte Carlo workflow with two-step `ctch_sinks` partitioning
- `randomly/`: advanced stochastic workflows using WWFR
- `theory/`: theory and diagnostics
- `data/`: example data and provenance scripts
