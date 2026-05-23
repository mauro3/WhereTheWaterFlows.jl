# Examples

Start Julia in the `examples/` folder with `julia --project` and run scripts via `include`.

To run all example scripts in one go:

```julia
include("run-all-examples.jl")
```

## Quick-start examples

- `wwf-simple.jl`: deterministic core routing with `WhereTheWaterFlows`
- `wwfs-simple.jl`: deterministic subglacial routing with `WhereTheWaterFlows.Subglacially`
- `wwfr-simple.jl`: stochastic routing using a Monte Carlo approach with `WhereTheWaterFlows.Randomly` using uncertainty in source term

## Topical examples

- `core/`: additional core WWF workflows
- `feedback_function/`: examples of custom `feedback_fn` models (melt and sediment transport)
- `subglacially/`: richer subglacial scenarios
- `randomly/`: advanced stochastic workflows
- `theory/`: theory and diagnostic scripts
- `data/`: input data and provenance scripts
