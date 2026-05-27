# Benchmark water routing

Benchmarking water-routing performance for [WhereTheWaterFlows.jl](https://github.com/mauro3/WhereTheWaterFlows.jl),
[Whitebox.jl](https://github.com/acgold/Whitebox.jl), [TopoToolbox](https://github.com/wschwanghart/topotoolbox), and [GRASS GIS](https://grass.osgeo.org/).

Benchmarks are to be taken with a grain of salt, as we may not run the other models in the best possible way. In particular the
multi-threaded runs do not show any/much gains, which maybe due to setup issues. PRs welcome!

## Latest generated benchmark table

[results.md](results.md)

Note: Whitebox writes intermediate files, and that I/O time is included in benchmark timings.

## Run scripts

### Setup

From this folder (`benchmarks/`):

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
```

### Run all laptop benchmarks

From this folder:

```bash
./run_all_benchmarks.sh
```

Useful environment variables:

- `JULIA_THREADS` (default `1`)
- `THREADS_LIST` (default `1,4,8`; comma-separated sweep for WWF/Whitebox)
- `WWF_RUNS` (default `6`, used for WWF small run)
- `WHITEBOX_RUNS` (default `6`)
- `WHITEBOX_LARGE_RUNS` (default `3`, used for Whitebox large run)
- `RUN_GRASS=auto|always|never` (default `auto`)
- `GRASS_BIN` (default `grass`)
- `GRASS_NPROCS` (default `1`, used only if your `r.watershed` supports `nprocs`)
- `GRASS_NPROCS_LIST` (default `$THREADS_LIST`; comma-separated sweep for GRASS)
- `GRASS_RUNS` (default `6`, used for GRASS small run)
- `GRASS_LARGE_RUNS` (default `3`, used for GRASS large run)
- `RUN_MATLAB=auto|always|never` (default `auto`)
- `MATLAB_BIN` (default `matlab`)
- `TOPO_THREADS` (default `1`)
- `TOPO_THREADS_LIST` (default `$THREADS_LIST`; comma-separated sweep for TopoToolbox)
- `TOPO_DATASET` (used by direct MATLAB calls; default `large`)
- `TOPO_DATASET_LIST` (default `small,large`; dataset sweep for TopoToolbox in `run_all_benchmarks.sh`)
- `TOPO_RUNS` (used by direct MATLAB calls; default `6`)
- `TOPO_SMALL_RUNS` (default `6`, used by `run_all_benchmarks.sh`)
- `TOPO_LARGE_RUNS` (default `3`, used by `run_all_benchmarks.sh`)
- `JULIA_BIN` (default `julia`)

Thread-sweep example for a single full pass:

```bash
THREADS_LIST=1,4,8 GRASS_NPROCS_LIST=1,4,8 TOPO_THREADS_LIST=1,4,8 TOPO_DATASET_LIST=small,large ./run_all_benchmarks.sh
```

Single-thread-count run example (2 threads only):

```bash
THREADS_LIST=2 GRASS_NPROCS_LIST=2 TOPO_THREADS_LIST=2 TOPO_DATASET_LIST=small,large ./run_all_benchmarks.sh
```

## WWF

Quick mock run:

```bash
julia --project run_wwf.jl --mode mock --dataset small --runs 1
```

Real run on the small DEM using 1 thread:

```bash
julia -t1 --project run_wwf.jl --mode real --dataset small --runs 6
```

Real run on the large DEM (default 1 thread, 3 runs):

```bash
julia --project run_wwf.jl --mode real --dataset large
```

## Whitebox

Whitebox is run through [Whitebox.jl](https://github.com/acgold/Whitebox.jl), which calls WhiteboxTools from Julia.
Each timed run includes DEM conditioning via `breach_depressions(fill_pits=false)` before D8 flow accumulation.

Quick mock run:

```bash
julia --project run_whitebox.jl --mode mock --dataset small --runs 1
```

Real run on the small DEM using 1 thread:

```bash
julia -t1 --project run_whitebox.jl --mode real --dataset small --runs 6
```

Real run on the large DEM (default 1 thread, 3 runs):

```bash
julia --project run_whitebox.jl --mode real --dataset large
```

Note: Whitebox large runs require `data_raw/swissalti3d_tilelarge.tif`.
Create it first via WWF, for example:

```bash
julia --project run_wwf.jl --mode real --dataset large --runs 1
```

## TopoToolbox

Install MATLAB and clone TopoToolbox:

```bash
git clone https://github.com/TopoToolbox/topotoolbox3.git topotoolbox3
```

Run from MATLAB prompt (after `cd` into `benchmarks/`):

```matlab
run('run_topotoolbox.m')
```

Run from terminal with defaults (`TOPO_THREADS=1`, `TOPO_DATASET=large`, `TOPO_RUNS=6`):

```bash
matlab -batch "run('run_topotoolbox.m')"
```

Run from terminal with overrides:

```bash
TOPO_THREADS=4 TOPO_DATASET=large TOPO_RUNS=3 matlab -batch "run('run_topotoolbox.m')"
```

Note: TopoToolbox large runs require `data_raw/swissalti3d_tilelarge.tif`.
Create it first via WWF, for example:

```bash
julia --project run_wwf.jl --mode real --dataset large --runs 1
```

## GRASS

GRASS is benchmarked through the command line using `r.watershed -s` (single-flow-direction / D8 mode).

Install options:

- macOS (Homebrew):

```bash
brew install grass
```

- Ubuntu/Debian:

```bash
sudo apt-get update
sudo apt-get install -y grass
```

Quick runs:

```bash
julia --project run_grass.jl --mode real --dataset small --runs 1
julia --project run_grass.jl --mode real --dataset large
```

Notes:

- Large GRASS runs require `data_raw/swissalti3d_tilelarge.tif` (create once via WWF large mode).
- To override binary or process count:

```bash
GRASS_BIN=grass GRASS_NPROCS=1 julia --project run_grass.jl --mode real --dataset small --runs 6
```
