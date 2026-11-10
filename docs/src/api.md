# [API Reference](@id APIReference)

Exported names are available after `using` the corresponding module. Names declared
`public` but not exported are supported through qualification, such as
`WWFS.GAMMA`. Internal functions are documented below for contributors, but their
interfaces are not supported public API.

High-level guides:

- [Tutorial](@ref)
- [Feedback Functionality](@ref FeedbackGuide)
- [Randomly](@ref RandomlyGuide)
- [Subglacially](@ref SubglaciallyGuide)

## Public API

### WhereTheWaterFlows

```@autodocs
Modules = [WhereTheWaterFlows]
Public = true
Private = false
Filter = f -> !(f isa Function && startswith(string(nameof(f)), "plt_"))
```

### Subglacially

```@autodocs
Modules = [WhereTheWaterFlows.Subglacially]
Public = true
Private = false
```

### Randomly

```@autodocs
Modules = [WhereTheWaterFlows.Randomly]
Public = true
Private = false
```

### Plotting (Makie extension)

Load Makie (or a backend such as CairoMakie) to enable these exported
`WhereTheWaterFlows.plt_*` functions. Their implementations and docstrings live
in the Makie extension, but the supported bindings belong to `WhereTheWaterFlows`.
Each recipe also has an exported mutating variant (`plt_dir!`, `plt_area!`,
`plt_sinks!`, `plt_bnds!`, `plt_catchments!`) for plotting into an existing axis.
`plt_it` has no mutating variant. Extension-local preprocessing helpers and
Makie-generated recipe types are internal.

```@docs
WhereTheWaterFlows.plt_dir
WhereTheWaterFlows.plt_area
WhereTheWaterFlows.plt_sinks
WhereTheWaterFlows.plt_bnds
WhereTheWaterFlows.plt_catchments
WhereTheWaterFlows.plt_it
```

## Internals

### WhereTheWaterFlows

```@autodocs
Modules = [WhereTheWaterFlows]
Public = false
Private = true
```

### Subglacially

```@autodocs
Modules = [WhereTheWaterFlows.Subglacially]
Public = false
Private = true
```

### Randomly

```@autodocs
Modules = [WhereTheWaterFlows.Randomly]
Public = false
Private = true
```
