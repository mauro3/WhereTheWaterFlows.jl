# Benchmark water routing

Benchmarking water-routing performance [WhereTheWaterFlows.jl](https://github.com/mauro3/WhereTheWaterFlows.jl),
[Whitebox.jl](https://github.com/acgold/Whitebox.jl), [TopoToolbox](https://github.com/wschwanghart/topotoolbox), and [GRASS GIS](https://grass.osgeo.org/).


The tests were performed on an Intel i7-8565U CPU using single-thread and 4-thread configurations
(i.e., 4 physical CPU cores) tested on two DEMs. The "small" DEM corresponds to a single swissALTI3D tile
(2000 x 2000: 4 million pixels, 3974 pits), while the "large" DEM corresponds to 4 swissALTI3D tile
(8000 × 8000: 64 million pixels, 11332 pits). All routing tested use a D8 algorithm.
Mean runtimes and standard deviations were computed from 6 repeated runs.

| Method | Software / functions | Threads | Runtime small DEM (s) | Runtime large DEM (s) |
| --- | --- | --- | --- | --- |
| **WWF** | Julia `waterflows` | 1 | **1.17 ± 0.02** | **20.31 ± 0.36** |
| **WWF** | Julia `waterflows` | 4 | **1.28 ± 0.02** | **22.27 ± 0.72** |
| **Whitebox** | Julia `d8_pointer` + `d8_flow_accumulation` | 1 | **2.88 ± 0.18** | **31.10 ± 0.77** |
| **Whitebox** | Julia `d8_pointer` + `d8_flow_accumulation` | 4 | **2.79 ± 0.22** | **29.89 ± 0.51** |
| **TopoToolbox** | MATLAB `FLOWobj`, `flowacc`, `drainagebasins` | 1 | **1.45 ± 0.06** | **33.20 ± 3.32** |
| **TopoToolbox** | MATLAB `FLOWobj`, `flowacc`, `drainagebasins` | 4 | **1.15 ± 0.02** | **33.05 ± 1.62** |
| **GRASS** | QGIS + GRASS `r.watershed` | unknown | **15.84** | **110.33** |


## Run scripts

### WWF

Clone this repository and activate the Julia environment:

```bash
git clone https://github.com/christopheogier/benchmark_water_routing.git
cd benchmark_water_routing
julia --project=scripts -e 'using Pkg; Pkg.instantiate()'
```

using 1 thread (in terminal at the repo location):
```bash
julia -t1 --project=scripts scripts/run_wwf.jl
```

using 4 thread:
```bash
julia -t4 --project=scripts scripts/run_wwf.jl
```

## Whitebox

Whitebox is run through [Whitebox.jl](https://github.com/acgold/Whitebox.jl), which calls WhiteboxTools from Julia.

Using 1 thread:

```bash
julia -t1 --project=scripts scripts/run_whitebox.jl
```

Using 1 thread:

```bash
julia -t4 --project=scripts scripts/run_whitebox.jl
```

### Topotoolbox

Install MATLAB and clone TopoToolbox with
```
git clone https://github.com/TopoToolbox/topotoolbox3.git
```

set manually `maxNumCompThreads(n)` with n=1 or 4 in `run_topotoolbox.m`, then:
```
run('scripts/run_topotoolbox.m')
```

### GRASS
Install QGIS.

The test on GRASS (`r.watershed`) was performed manually through the QGIS graphical interface and is therefore not directly reproducible using the scripts provided in this repository.

