# [Feedback Functionality (advanced usage)](@id FeedbackGuide)

`waterflows` can route one or several *extensive* quantities downstream. By
providing `feedback_fn`, you can update those quantities at each cell after all
upstream contributions have been accumulated but before the total is passed to
the downstream receiver.

This is useful for coupled process models where routed water modifies itself or
other transported quantities, for example:

- dissipation melt generation in a subglacial setting (as is done in WWFS),
- sediment transport capacity limits,
- simple reactive tracer/source-sink terms.

## When do you need `feedback_fn`?

There are three distinct modes of WWF use:

| Mode | When to use it |
|---|---|
| **Pure routing** | Omit `feedback_fn`; supply source quantities via `cellarea`. |
| **Local source injection** | The source at each cell is computed from a non-routed auxiliary field (e.g. snow depth, rainfall intensity). Close over the auxiliary array and inject into `uparea` inside the callback; the auxiliary field is never accumulated itself. |
| **State-dependent modification** | The quantity leaving a cell depends on how much has arrived from upstream (capacity limits, storage release, reactive loss). Only a callback — not a pre-computed raster — can express this, because the rule must be evaluated at traversal time. |

Real process models often combine all three.

### Local source injection: snow melt example

A linear snow-melt model where discharge is driven by local melt (proportional
to snow depth) rather than a uniform source illustrates the second mode:

```julia
snowdepth = ...      # auxiliary field — read inside callback, never routed
k         = 0.1      # melt rate constant

function snowmelt_feedback(uparea_Q, ij, dir)
    local_melt = k * snowdepth[ij]   # local-only computation
    return uparea_Q + local_melt     # inject into routed discharge
end

out = waterflows(dem, zeros(size(dem)); feedback_fn=snowmelt_feedback)
```

`cellarea` is zero — no source is pre-specified.  The entire water input is
computed locally inside the callback from the non-routed `snowdepth` field.
If the melt rate also depended on incoming discharge (e.g. warm-water melting),
that would be the third mode, and `uparea_Q` would appear in the melt formula.

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
once all upstream cells have been visited, in upstream-to-downstream order.

## Minimal example

```@example feedback
using WhereTheWaterFlows
using Statistics

n = 80
x = range(-pi, pi, length=n)
dem = sin.(x) .* cos.(x')

# route uniform source and clip if uparea goes above 1000
clip_more_than_1000(uparea, ij, dir) = min(uparea, 1000)

out = waterflows(dem, ones(size(dem)); feedback_fn=clip_more_than_1000)
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

Note that this is essentially how WWFS calculates extra subglacial melt and the
code of WWFS can serve as another example. 

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

See also: [Examples](@ref ExamplesPage) and [API Reference](@ref APIReference).
