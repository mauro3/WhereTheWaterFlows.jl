# Examples

Start Julia in the `examples/` folder with `julia --project` and run scripts via `include`.

To run all example scripts in one go:

```julia
include("run-all-examples.jl")
```

## Quick-start examples

- `wwf-simple.jl`: deterministic core routing with `WhereTheWaterFlows`
- `wwfs-simple.jl`: deterministic subglacial routing with `WhereTheWaterFlows.Subglacially`
- `wwfr-simple.jl`: stochastic routing with `WhereTheWaterFlows.Randomly` using uncertainty in source term

These three are used as smoke examples in the package tests.

## Topical examples

- `core/`: additional core WWF workflows
- `feedback_function/`: examples of custom `feedback_fn` models (melt and sediment transport)
- `subglacially/`: richer subglacial scenarios
- `randomly/`: advanced stochastic workflows
- `theory/`: theory and diagnostic scripts
- `data/`: input data and provenance scripts

## Suggested learning path

1. `wwf-simple.jl`
2. `wwfs-simple.jl`
3. `wwfr-simple.jl`
4. then pick a topical folder (`feedback_function/`, `subglacially/`, `randomly/`)
