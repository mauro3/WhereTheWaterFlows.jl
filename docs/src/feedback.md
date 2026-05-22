# [Feedback Functionality](@id FeedbackGuide)

`waterflows` can route one or several *extensive* quantities downstream. By
providing `feedback_fn`, you can update those quantities at each cell before
they are passed to the downstream receiver.

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

- `uparea`: the accumulated routed quantity/quantities at `ij`
- `ij`: current `CartesianIndex`
- `dir`: flow-direction matrix
- return value must match the shape used by `cellarea`

If `cellarea` is a single array, `uparea` is a scalar. If `cellarea` is a
tuple of arrays, `uparea` is a tuple with matching structure.

In all cases, `feedback_fn` is called after upstream accumulation at the cell
and before routing to the downstream receiver.

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
    # local update (placeholder): damp tracer with discharge
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
- For slope-based laws, handle sink cells (`dir[ij] >= SINK`) explicitly.
- Start simple, then add complexity once baseline behavior is validated.

See also: [Examples](@ref) and [API Reference](@ref).
