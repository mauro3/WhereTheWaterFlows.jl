---
title: 'WhereTheWaterFlows.jl: Hydrological flow routing on digital elevation models in Julia'
tags:
  - Julia
  - hydrology
  - glaciology
  - digital elevation model
  - flow routing
  - D8
  - uncertainty quantification
authors:
  - name: Mauro A. Werder
    orcid: 0000-0003-0137-9377
    corresponding: true
    affiliation: "1, 2"
  - name: Christophe Ogier
    orcid: 0000-0002-5526-6071
    affiliation: "1, 2"
affiliations:
  - name: Laboratory of Hydraulics, Hydrology and Glaciology (VAW), ETH Zurich, Zurich, Switzerland
    index: 1
    ror: "05a28rw58"
  - name: Swiss Federal Institute for Forest, Snow and Landscape Research (WSL), bâtiment ALPOLE, Sion, Switzerland
    index: 2
    ror: "04bs5yc70"
date: 28 May 2026
bibliography: paper.bib
---

# Summary

`WhereTheWaterFlows.jl` (`WWF`) is a Julia package for computing water flow paths and
drainage basins on digital elevation models (DEMs) or other hydropotential surfaces.
Given a DEM, the package determines *where the water flows* across the
landscape, accumulating upslope area and delineating catchment basins.
Besides this classic, deterministic routing on surface topographies, WWF also includes routines
to route water in subglacial settings as well as uncertainty quantification tools based on a
Monte Carlo approach.

Key outputs include the flow-direction field, the
upslope contributing area or discharge, the
depression-filled DEM, a labelled catchment-basin map, and the location of all sinks.
Plotting functions for result visualisation are included. The package is
registered in the Julia General registry, has a full test suite and complete
[online documentation](https://mauro3.github.io/WhereTheWaterFlows.jl/stable).

WWF implements the widely-used D8 single-flow-direction algorithm, which routes water from each grid cell to its
lowest-elevation neighbour among the eight surrounding cells, together with a breach-type
depression-filling algorithm to handle closed depressions that would otherwise
trap flow [@OCallaghanExtractionDrainageNetworks1984]. The underlying graph traversal follows the O(n) recursive strategy
of @BraunVeryEfficientImplicit2013, giving linear scaling with the number of grid cells.
WWF performs comparably to or better than existing routing software, as shown in
our [benchmarks](https://github.com/mauro3/WhereTheWaterFlows.jl/blob/master/benchmarks/README.md).

# Statement of need

Hydrological flow routing on DEMs is a fundamental operation in geosciences:
it underpins catchment delineation, runoff modelling, subglacial hydrology,
and landscape-evolution studies. Software implementations of flow routing date
back to @OCallaghanExtractionDrainageNetworks1984 and, nowadays, many implementations exist. However,
`WWF` is the first, and to our knowledge, only native Julia package for this task.

WWF is suitable for any domain where DEM-based flow analysis is required.
However, the package has been designed with glaciological applications in mind where it
has been used in six publications to date, see section [Research impact](#research-impact).
The glaciological focus required the implementation of glacier specific aspects, such as Shreve
potential routing [@ShreveMovementWaterGlaciers1972] in the `Subglacially` submodule, as well as the development
of uncertainty quantification via Monte Carlo methods (the `Randomly` submodule). The latter
aspect is needed in glaciology as uncertainties in glacier bed DEMs and Shreve-potential are
large and thus warrant quantification.

# State of the field

Several mature tools exist for hydrological flow routing on DEMs, examples include:

- **GRASS GIS** (command line/C/Python/QGIS), the component `r.watershed` implements
  several water routing and related algorithms. [@GRASSDevelopmentTeamGRASS2026]
- **TauDEM** (command line/ArcGIS) provides command-line tools designed for large-scale catchment
  analysis. [@TesfaExtractionHydrologicalProximity2011]
- **TopoToolbox** (MATLAB/Python) provides widely used functionalities
  for topographic analysis including flow routing. [@SchwanghartShortCommunicationTopoToolbox2014]
- **Whitebox Tools** (Rust/Python) provides an extensive geomorphometry
  toolkit as a Rust library with Python bindings and unofficial [Julia bindings](https://github.com/acgold/Whitebox.jl).
  [@LindsayWhiteboxGATCase2016] 
- **RichDEM** (Python/C++) offers a broad set of flow-routing methods. [@BarnesPriorityfloodOptimalDepressionfilling2014]

To our knowledge, none of the existing routing software packages implement specialised
subglacial routing functionality or have uncertainty quantification capabilities.

# Software design

The central function is `waterflows(dem)`, which accepts a two-dimensional
array of elevation values and returns a named tuple of outputs such as upstream
area and catchments.
Boundary conditions (open or closed edges, fixed sinks) and other options can be specified via
keyword arguments.
A novel feature of WWF is
the feedback-function (via the `feedback_fn` argument of `waterflows`) which is called on each cell before
routing its water downstream, enabling feedbacks to occur. This allows for
extensive customisations, as demonstrated in the provided examples (see `examples/` directory)
on sediment transport and subglacial melt calculations.

For subglacial routing, the main function is `waterflows_subglacial(surfdem, beddem, dx, f)`
with surface and bed DEMs, grid spacing and water pressure as flotation fraction as inputs.
In the subglacial setting, the `beddem` and `f` input fields have large uncertainties associated with them.
This spurred the development of the Monte Carlo based uncertainty quantification in WWF. The main function
there is `map_mc(model, sample, reduce!, n)` which runs the `model` function (i.e. the water routing)
on `n` samples of surface, bed and flotation fractions fields and reduces the results with `reduce!`
(to avoid huge data volumes). The uncertainty in the input fields is modelled using Gaussian Random Fields,
which provide spatially correlated random noise. The random fields are generated using a FFT-based method [@RassEfficientParallelRandom2019]. However, the software design allows to couple
more sopisticated, external geostatistical models to generete the uncertainty fields.


# Example

The code and figure below show an example of a drainage network calculation on a small synthetic DEM:
```julia
using WhereTheWaterFlows, CairoMakie
n   = 200
x = y = range(-π, π, length=n)
dem = sin.(x) .* cos.(y') .+ 0.05 * rand(n, n)

out = waterflows(dem)
plt_area(x, y, out.area)
```
![Upslope contributing area on a synthetic DEM from running the provided code example. Brighter colours indicate
larger contributing areas, tracing the drainage network emerging from the
synthetic topography.\label{fig:area}](upslope_area.png)

# Research impact

`WhereTheWaterFlows.jl` has been used in published glaciological research:
@MalczykConstraintsSubglacialMelt2023 used it to route subglacial melt water and constrain melt
  fluxes beneath an ice sheet from observations of active subglacial lake
  recharge. The study employed WWF's uncertainty quantification features.
@DelaneyModelingSpatiallyDistributed2023 built a 2D subglacial sediment transport model on top of WWF.
@OgierDefinitionFormationRupture2025 employed flow routing to analyse the likelihood of subglacial water pockets in Alpine glaciers, contributing to a new
  inventory for the Swiss Alps.
@HorganWestAntarcticGroundingzone2025 used flow routing to delineate the subglacial catchment of the Kamb Ice Stream.
@WashamOceanicVolcanicHeat2026 demonstrated that volcanic heat sources within the Kamb catchment likely
contribute to enhanced subglacial melt.
@OgierPotentialGlacierContributions2026 (in review) used the package to assess the potential contribution of subglacial water reservoirs to the catastrophic 2024 La Bérarde flood.

WWF is used for teaching in the lecture courses "Physics of Glaciers" at ETH Zurich (Switzerland) and
"Introduction to geoscientific programming" at Uni Mainz (Germany).
The software is archived on Zenodo [@MauroA.WhereTheWaterFlowsjl2024] and has
accumulated 20+ GitHub stars and 5+ forks since its public release in 2019.
It has four contributors and 200+ commits across 20+ tagged releases.

# AI usage disclosure

AI tools from Anthropic and OpenAI are used in the development of WWF. Most of
the code predates 2022, the advent of capable LLMs, and was thus written by hand.
This also applies to the `Randomly` and `Subglacially` submodules which were previously
developed in a separate, private repository
and only now merged into WWF. Recently, code refactoring was helped by LLMs.
Much of the documentation, examples and benchmarks were generated
with the help of LLMs. The tools employed to date are Claude Sonnet 4.6
and GPT-5.3-Codex. All AI generated content was reviewed by a human.

The first draft of this manuscript was created by Claude Sonnet 4.6.
The tool was used for initial text generation
and structural organisation of the paper sections. The text was then reviewed,
heavily edited and rewritten by the two authors.

# Acknowledgements

Development of WWF was partly funded by the European Space Agency's project 4DAntarctica (ESA: Grant 4000128611/19/I-DT)
and by the Swiss National Science Foundation's project DIWING (grant nr. 212061).

# References
