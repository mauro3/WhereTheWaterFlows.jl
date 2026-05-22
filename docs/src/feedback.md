# [Feedback Functionality](@id FeedbackGuide)

`waterflows` can route one or several *extensive* quantities downstream. By
providing `feedback_fn`, you can update those quantities at each cell after all
upstream contributions have been accumulated but before the total is passed to
the downstream receiver.

This is useful for coupled process models where routed water modifies itself or
other transported quantities, for example:

- dissipation melt generation,
- sediment transport capacity limits,
- simple reactive tracer/source-sink terms.

## Feedback contract

The callback signature is:

```julia
feedback_fn(uparea, ij, dir) -> new_uparea
```

- `uparea`: accumulated quantity/quantities arriving at `ij` from all upstream cells, plus `cellarea[ij]` itself
- `ij`: current `CartesianIndex`
- `dir`: the full flow-direction matrix (read-only)
- return value: the updated quantity to be passed downstream; must have the same type and shape as `uparea`

When `cellarea` is a single array, `uparea` is a scalar.  When `cellarea` is a
tuple of arrays, `uparea` is a tuple of scalars with matching structure.

The callback is invoked at every active cell (i.e. cells that are not barriers)
once all upstream cells have been visited, in up-to-downstream order.

## Minimal scalar feedback

```@example feedback
using WhereTheWaterFlows
using Statistics

n = 80
x = range(-pi, pi, length=n)
dem = sin.(x) .* cos.(x')

# route uniform source and clip negative values (toy example)
clip_nonnegative(uparea, ij, dir) = max(uparea, 0.0)

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

## Typical implementation skeleton

```@example feedback
using WhereTheWaterFlows

n = 40
x = range(-pi, pi, length=n)
dem = sin.(x) .* cos.(x')

cellarea = (ones(size(dem)), fill(0.1, size(dem)))

function toy_feedback(uparea, ij, dir)
    Q, tracer = uparea
    d = dir[ij]
    if d >= WhereTheWaterFlows.SINK
        return (Q, tracer)
    end
    # damp tracer with accumulated discharge (placeholder)
    tracer_new = tracer / (1 + 0.01 * Q)
    return (Q, tracer_new)
end

out = waterflows(dem, cellarea; feedback_fn=toy_feedback)
mean(out.area[2])
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

## Common patterns

- **Clipping or limiting**: enforce positivity or capacity constraints.
- **Production terms**: add melt/erosion/reaction source terms.
- **Diagnostic fields**: route extra arrays for bookkeeping outputs.
- **Coupled extensive fields**: route `(water, tracer)` or `(water, sediment)`.

## Practical tips

- Use only extensive/additive routed quantities.
- Keep units consistent with your `cellarea` definition.
- For slope-based laws, check `dir[ij] >= WhereTheWaterFlows.SINK` to skip
  cells with no valid downstream direction.
- Check the return value at sink/barrier cells explicitly if the process law
  would produce nonsensical results there.
- Start simple, then add complexity once baseline behaviour is validated.

See also: [Examples](@ref) and [API Reference](@ref).
