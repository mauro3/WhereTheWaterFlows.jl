# Benchmark Results

Generated: 2026-05-28 10:47:11
Platform: Linux
CPU: AMD Ryzen AI 7 PRO 350 w/ Radeon 860M
Logical threads: 16

## WWF

| method | dataset | dem_size | features | nthreads | nruns | runtime_mean_s | runtime_std_s |
| --- | --- | --- | --- | --- | --- | --- | --- |
| wwf | small | 2000x2000 | waterflows (real, small) | 1 | 6 | 0.491 | 0.008 |
| wwf | small | 2000x2000 | waterflows (real, small) | 4 | 6 | 0.488 | 0.008 |
| wwf | small | 2000x2000 | waterflows (real, small) | 8 | 6 | 0.503 | 0.010 |
| wwf | large | 8000x8000 | waterflows (real, large) | 1 | 3 | 7.733 | 0.080 |
| wwf | large | 8000x8000 | waterflows (real, large) | 4 | 3 | 8.043 | 0.087 |
| wwf | large | 8000x8000 | waterflows (real, large) | 8 | 3 | 7.879 | 0.002 |

## GRASS

| method | dataset | dem_size | features | nthreads | nruns | runtime_mean_s | runtime_std_s |
| --- | --- | --- | --- | --- | --- | --- | --- |
| grass | small | 2000x2000 | r.watershed -s (D8, real, small, nprocs=unsupported) | missing | 6 | 1.203 | 0.008 |
| grass | large | 8000x8000 | r.watershed -s (D8, real, large, nprocs=unsupported) | missing | 3 | 18.513 | 0.064 |

## Whitebox

| method | dataset | dem_size | features | nthreads | nruns | runtime_mean_s | runtime_std_s |
| --- | --- | --- | --- | --- | --- | --- | --- |
| whitebox | small | 2000x2000 | breach_depressions(fill_pits=false) + d8_flow_accumulation (real, small) | 1 | 6 | 2.489 | 0.011 |
| whitebox | small | 2000x2000 | breach_depressions(fill_pits=false) + d8_flow_accumulation (real, small) | 4 | 6 | 2.497 | 0.006 |
| whitebox | small | 2000x2000 | breach_depressions(fill_pits=false) + d8_flow_accumulation (real, small) | 8 | 6 | 2.488 | 0.023 |
| whitebox | large | 8000x8000 | breach_depressions(fill_pits=false) + d8_flow_accumulation (real, large) | 1 | 3 | 41.378 | 0.098 |
| whitebox | large | 8000x8000 | breach_depressions(fill_pits=false) + d8_flow_accumulation (real, large) | 4 | 3 | 41.431 | 0.225 |
| whitebox | large | 8000x8000 | breach_depressions(fill_pits=false) + d8_flow_accumulation (real, large) | 8 | 3 | 41.591 | 0.102 |

## TopoToolbox

| method | dataset | dem_size | features | nthreads | nruns | runtime_mean_s | runtime_std_s |
| --- | --- | --- | --- | --- | --- | --- | --- |
| topotoolbox | small | 2000x2000 | FLOWobj + flowacc + drainagebasins | 1 | 6 | 0.760 | 0.020 |
| topotoolbox | small | 2000x2000 | FLOWobj + flowacc + drainagebasins | 4 | 6 | 0.533 | 0.022 |
| topotoolbox | small | 2000x2000 | FLOWobj + flowacc + drainagebasins | 8 | 6 | 0.517 | 0.027 |
| topotoolbox | large | 8000x8000 | FLOWobj + flowacc + drainagebasins | 1 | 3 | 15.836 | 0.096 |
| topotoolbox | large | 8000x8000 | FLOWobj + flowacc + drainagebasins | 4 | 3 | 11.979 | 0.087 |
| topotoolbox | large | 8000x8000 | FLOWobj + flowacc + drainagebasins | 8 | 3 | 11.703 | 0.112 |
