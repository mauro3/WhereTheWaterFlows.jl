# [Feedback Functionality](@id FeedbackGuide)

`waterflows` can route one or several *extensive* quantities downstream. By
providing a `feedback_fn`, you can update those quantities at each cell after all
upstream contributions have been accumulated but before the total is passed to
the downstream receiver.

This is useful for coupled process models where routed water modifies itself or
other transported quantities, for example:

- dissipation melt generation in a subglacial setting,
- sediment transport,
- simple reactive tracer/source-sink terms.
- local water storage (e.g. a linear storage model at each cell)

Potentially, other process models, for instance a glacier melt model, could be implemented as
`feedback_fn` opening a way to extend WWF with other models.

## Feedback contract

The callback signature is:

```julia
feedback_fn(uparea_incoming, ij, dir) -> uparea_to_be_routed
```

- `uparea_incoming`: accumulated quantity/quantities arriving at `ij` from all upstream cells, plus `cellarea[ij]` itself
- `ij`: current `CartesianIndex`
- `dir`: the full flow-direction matrix (read-only)
- return value: the updated quantity to be passed downstream; must have the same type and shape as `uparea_incoming`

When `cellarea` is a single array, `uparea` is a scalar.  When `cellarea` is a
tuple of arrays, `uparea` is a tuple of scalars with matching structure.

The callback is invoked at every active cell (i.e. cells that are not barriers)
once all upstream cells have been visited, in up-to-downstream order.

## Minimal example

```@example feedback
using WhereTheWaterFlows
using Statistics

n = 80
x = range(-pi, pi, length=n)
dem = sin.(x) .* cos.(x')

# route uniform source and clip if uparea goes above 1000
clip_more_than_1000(uparea, ij, dir) = min(uparea, 1000)

out = waterflows(dem, ones(size(dem)); feedback_fn=clip_nonnegative)
maximum(out.area)
```

## Multi-field feedback (water + sediment)

The `examples/feedback_function/sediment-transport-mpm.jl` script shows a
simple MPM-like sediment capacity closure using two routed fields:

- water discharge `Q`
- sediment discharge `Qs`

with a stream-width scaling `W ~ Q^b`.

Conceptually:

```julia
Qs_cap = capacity(Q, local_slope)
Qs_out = min(Qs_in, Qs_cap)
```

and the feedback returns `(Q, Qs_out)`.

You can run that full example directly from the examples environment:

```julia
include("feedback_function/sediment-transport-mpm.jl")
```

## Multi-field feedback (subglacial melt)

`examples/feedback_function/subglacial-routing-feedback.jl` illustrates a
three-field setup that routes:

- total discharge,
- extra discharge due to melt,
- local melt rate diagnostic.

This pattern is especially useful when you need one or more additional fields
for diagnostics while keeping the physically relevant transport state explicit.

Run the full example with:

```julia
include("feedback_function/subglacial-routing-feedback.jl")
```

Note that this is essentially how WWFS calculates extra subglacial melt. 

## Tips and tricks

Common patterns

- **Clipping or limiting**: enforce positivity or capacity constraints.
- **Production terms**: add melt/erosion/reaction source terms.
- **Diagnostic fields**: route extra arrays for bookkeeping outputs.
- **Coupled extensive fields**: route `(water, tracer)` or `(water, sediment)`.

Tips

- Use only extensive/additive routed quantities.
- Keep units consistent with your `cellarea` definition.
- For slope-based laws, check `dir[ij] >= WhereTheWaterFlows.SINK` to skip
  cells with no valid downstream direction.
- Check the return value at sink/barrier cells explicitly if the process law
  would produce nonsensical results there.

See also: [Examples](@ref) and [API Reference](@ref).
