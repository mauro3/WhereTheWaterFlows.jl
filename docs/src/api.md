# API Reference

## Public API

```@docs
waterflows
fill_dem
```

High-level guides:

- [Tutorial](@ref)
- [Feedback Functionality](@ref FeedbackGuide)
- [Randomly](@ref RandomlyGuide)
- [Subglacially](@ref SubglaciallyGuide)

## Post-processing

These functions are exported but are called separately after `waterflows`:

```@docs
catchment
catchments
catchment_flux
prune_catchments
```

## Internals

These functions are not part of the public API but are documented for
contributors and users who want to customise the routing pipeline.

```@docs
WhereTheWaterFlows.d8dir_feature
WhereTheWaterFlows.flowrouting_catchments
WhereTheWaterFlows.drainpits!
WhereTheWaterFlows.make_boundaries
WhereTheWaterFlows.make_flowfeatures
WhereTheWaterFlows.dir2ind
WhereTheWaterFlows.dir2vec
WhereTheWaterFlows.ind2dir
WhereTheWaterFlows.diagonal_fac
WhereTheWaterFlows.flowsinto
WhereTheWaterFlows.on_outer_boundary
WhereTheWaterFlows.iterate_D9
WhereTheWaterFlows.showme
WhereTheWaterFlows.dirnums
WhereTheWaterFlows.cartesian
WhereTheWaterFlows._flow_from_to!
```

## Subglacially

```@autodocs
Modules = [WhereTheWaterFlows.Subglacially]
```

## Randomly

```@autodocs
Modules = [WhereTheWaterFlows.Randomly]
```
