#!/usr/bin/env sh
# Runs all local laptop benchmarks.
set -eu
set -o pipefail 2>/dev/null || true

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
BENCH_DIR="$ROOT_DIR/benchmarks"
OUT_DIR="$BENCH_DIR/outputs"

JULIA_BIN="${JULIA_BIN:-julia}"
JULIA_THREADS="${JULIA_THREADS:-1}"
THREADS_LIST="${THREADS_LIST:-1,4,8}"
WWF_RUNS="${WWF_RUNS:-6}"
WHITEBOX_RUNS="${WHITEBOX_RUNS:-6}"
WHITEBOX_LARGE_RUNS="${WHITEBOX_LARGE_RUNS:-3}"
GRASS_RUNS="${GRASS_RUNS:-6}"
GRASS_LARGE_RUNS="${GRASS_LARGE_RUNS:-3}"
GRASS_BIN="${GRASS_BIN:-grass}"
GRASS_NPROCS="${GRASS_NPROCS:-1}"
GRASS_NPROCS_LIST="${GRASS_NPROCS_LIST:-$THREADS_LIST}"
RUN_GRASS="${RUN_GRASS:-auto}" # auto|always|never

RUN_MATLAB="${RUN_MATLAB:-auto}" # auto|always|never
MATLAB_BIN="${MATLAB_BIN:-matlab}"
TOPO_THREADS="${TOPO_THREADS:-1}"
TOPO_THREADS_LIST="${TOPO_THREADS_LIST:-$THREADS_LIST}"
TOPO_DATASET="${TOPO_DATASET:-large}"
TOPO_DATASET_LIST="${TOPO_DATASET_LIST:-small,large}"
TOPO_RUNS="${TOPO_RUNS:-6}"
TOPO_SMALL_RUNS="${TOPO_SMALL_RUNS:-6}"
TOPO_LARGE_RUNS="${TOPO_LARGE_RUNS:-3}"

LARGE_DEM_FILE="$BENCH_DIR/data_raw/swissalti3d_tilelarge.tif"

cleanup() {
  if [ -f "$LARGE_DEM_FILE" ]; then
    echo "==> Cleaning up generated large DEM"
    rm -f "$LARGE_DEM_FILE"
  fi
}

trap cleanup EXIT

validate_positive_int() {
  case "$1" in
    ''|*[!0-9]*)
      return 1
      ;;
    0)
      return 1
      ;;
    *)
      return 0
      ;;
  esac
}

for th in $(printf '%s' "$THREADS_LIST" | tr ',' ' '); do
  validate_positive_int "$th" || {
    echo "ERROR: THREADS_LIST must be comma-separated positive integers" >&2
    exit 1
  }
done

for th in $(printf '%s' "$GRASS_NPROCS_LIST" | tr ',' ' '); do
  validate_positive_int "$th" || {
    echo "ERROR: GRASS_NPROCS_LIST must be comma-separated positive integers" >&2
    exit 1
  }
done

for th in $(printf '%s' "$TOPO_THREADS_LIST" | tr ',' ' '); do
  validate_positive_int "$th" || {
    echo "ERROR: TOPO_THREADS_LIST must be comma-separated positive integers" >&2
    exit 1
  }
done

for ds in $(printf '%s' "$TOPO_DATASET_LIST" | tr ',' ' '); do
  case "$ds" in
    small|large)
      ;;
    *)
      echo "ERROR: TOPO_DATASET_LIST must contain only 'small' and/or 'large'" >&2
      exit 1
      ;;
  esac
done

echo "==> Instantiating benchmarks environment"
"$JULIA_BIN" --project="$BENCH_DIR" -e 'using Pkg; Pkg.instantiate()'

echo "==> Removing previous benchmark CSV outputs"
rm -f "$OUT_DIR/wwf/benchmark_wwf.csv"
rm -f "$OUT_DIR/whitebox/benchmark_whitebox.csv"
rm -f "$OUT_DIR/topotoolbox/benchmark_topotoolbox.csv"
rm -f "$OUT_DIR/grass/benchmark_grass.csv"

for julia_threads in $(printf '%s' "$THREADS_LIST" | tr ',' ' '); do
  echo "==> Running WWF benchmark (small DEM, threads=$julia_threads)"
  JULIA_NUM_THREADS="$julia_threads" "$JULIA_BIN" --project="$BENCH_DIR" "$BENCH_DIR/run_wwf.jl" --mode real --dataset small --runs "$WWF_RUNS"

  echo "==> Running WWF benchmark (large DEM, threads=$julia_threads)"
  JULIA_NUM_THREADS="$julia_threads" "$JULIA_BIN" --project="$BENCH_DIR" "$BENCH_DIR/run_wwf.jl" --mode real --dataset large

  echo "==> Running Whitebox benchmark (small DEM, threads=$julia_threads)"
  JULIA_NUM_THREADS="$julia_threads" "$JULIA_BIN" --project="$BENCH_DIR" "$BENCH_DIR/run_whitebox.jl" --mode real --dataset small --runs "$WHITEBOX_RUNS"

  echo "==> Running Whitebox benchmark (large DEM, threads=$julia_threads)"
  JULIA_NUM_THREADS="$julia_threads" "$JULIA_BIN" --project="$BENCH_DIR" "$BENCH_DIR/run_whitebox.jl" --mode real --dataset large --runs "$WHITEBOX_LARGE_RUNS"
done

run_grass_benchmark() {
  grass_nprocs="$1"

  if ! command -v "$GRASS_BIN" >/dev/null 2>&1; then
    return 1
  fi

  echo "==> Running GRASS benchmark (small DEM, nprocs=$grass_nprocs)"
  GRASS_BIN="$GRASS_BIN" GRASS_NPROCS="$grass_nprocs" "$JULIA_BIN" --project="$BENCH_DIR" "$BENCH_DIR/run_grass.jl" --mode real --dataset small --runs "$GRASS_RUNS"

  echo "==> Running GRASS benchmark (large DEM, nprocs=$grass_nprocs)"
  GRASS_BIN="$GRASS_BIN" GRASS_NPROCS="$grass_nprocs" "$JULIA_BIN" --project="$BENCH_DIR" "$BENCH_DIR/run_grass.jl" --mode real --dataset large --runs "$GRASS_LARGE_RUNS"

  return 0
}

grass_supports_nprocs() {
  if ! command -v "$GRASS_BIN" >/dev/null 2>&1; then
    return 2
  fi

  probe_dir="$(mktemp -d /tmp/wwf_grass_probe.XXXXXX)"
  probe_loc="$probe_dir/location"
  probe_help="$($GRASS_BIN -c EPSG:4326 "$probe_loc" --exec r.watershed --help 2>&1 || true)"
  rm -rf "$probe_dir"

  case "$probe_help" in
    *nprocs*)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

grass_supports_nprocs_state="unknown"
if command -v "$GRASS_BIN" >/dev/null 2>&1; then
  if grass_supports_nprocs; then
    grass_supports_nprocs_state="yes"
  else
    grass_supports_nprocs_state="no"
  fi
fi

case "$RUN_GRASS" in
  always)
    if [ "$grass_supports_nprocs_state" = "no" ]; then
      set -- $(printf '%s' "$GRASS_NPROCS_LIST" | tr ',' ' ')
      first_grass_nprocs="$1"
      echo "==> GRASS r.watershed does not support nprocs; running once with nprocs=$first_grass_nprocs"
      run_grass_benchmark "$first_grass_nprocs" || {
        echo "ERROR: RUN_GRASS=always but GRASS is not available." >&2
        exit 1
      }
    else
      for grass_nprocs in $(printf '%s' "$GRASS_NPROCS_LIST" | tr ',' ' '); do
        run_grass_benchmark "$grass_nprocs" || {
          echo "ERROR: RUN_GRASS=always but GRASS is not available." >&2
          exit 1
        }
      done
    fi
    ;;
  never)
    echo "==> Skipping GRASS benchmark (RUN_GRASS=never)"
    ;;
  auto)
    if [ "$grass_supports_nprocs_state" = "no" ]; then
      set -- $(printf '%s' "$GRASS_NPROCS_LIST" | tr ',' ' ')
      first_grass_nprocs="$1"
      if run_grass_benchmark "$first_grass_nprocs"; then
        echo "==> GRASS r.watershed does not support nprocs; skipped remaining GRASS thread sweep values"
      else
        echo "==> Skipping GRASS benchmark (GRASS not found)"
      fi
    else
      for grass_nprocs in $(printf '%s' "$GRASS_NPROCS_LIST" | tr ',' ' '); do
        if run_grass_benchmark "$grass_nprocs"; then
          :
        else
          echo "==> Skipping GRASS benchmark (GRASS not found)"
          break
        fi
      done
    fi
    ;;
  *)
    echo "ERROR: RUN_GRASS must be one of: auto, always, never" >&2
    exit 1
    ;;
esac

run_matlab_benchmark() {
  topo_threads="$1"
  topo_dataset="$2"
  topo_runs="$3"
  topotoolbox_dir="$BENCH_DIR/topotoolbox3"

  if ! command -v "$MATLAB_BIN" >/dev/null 2>&1; then
    return 1
  fi

  if [ ! -d "$topotoolbox_dir" ]; then
    return 2
  fi

  echo "==> Running TopoToolbox benchmark (MATLAB, dataset=$topo_dataset, threads=$topo_threads, runs=$topo_runs)"
  TOPO_THREADS="$topo_threads" TOPO_DATASET="$topo_dataset" TOPO_RUNS="$topo_runs" "$MATLAB_BIN" -batch "run('$BENCH_DIR/run_topotoolbox.m')"
  return 0
}

case "$RUN_MATLAB" in
  always)
    for topo_dataset in $(printf '%s' "$TOPO_DATASET_LIST" | tr ',' ' '); do
      case "$topo_dataset" in
        small)
          topo_runs="$TOPO_SMALL_RUNS"
          ;;
        large)
          topo_runs="$TOPO_LARGE_RUNS"
          ;;
      esac
      for topo_threads in $(printf '%s' "$TOPO_THREADS_LIST" | tr ',' ' '); do
        run_matlab_benchmark "$topo_threads" "$topo_dataset" "$topo_runs" || {
          echo "ERROR: RUN_MATLAB=always but MATLAB/TopoToolbox setup is missing." >&2
          exit 1
        }
      done
    done
    ;;
  never)
    echo "==> Skipping TopoToolbox benchmark (RUN_MATLAB=never)"
    ;;
  auto)
    for topo_dataset in $(printf '%s' "$TOPO_DATASET_LIST" | tr ',' ' '); do
      case "$topo_dataset" in
        small)
          topo_runs="$TOPO_SMALL_RUNS"
          ;;
        large)
          topo_runs="$TOPO_LARGE_RUNS"
          ;;
      esac
      for topo_threads in $(printf '%s' "$TOPO_THREADS_LIST" | tr ',' ' '); do
        if run_matlab_benchmark "$topo_threads" "$topo_dataset" "$topo_runs"; then
          :
        else
          echo "==> Skipping TopoToolbox benchmark (MATLAB or benchmarks/topotoolbox3 not found)"
          break 2
        fi
      done
    done
    ;;
  *)
    echo "ERROR: RUN_MATLAB must be one of: auto, always, never" >&2
    exit 1
    ;;
esac

echo "==> Writing benchmark markdown report"
"$JULIA_BIN" --project="$BENCH_DIR" "$BENCH_DIR/write_results_md.jl"

echo "==> Done."
