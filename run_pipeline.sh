#!/usr/bin/env bash
# run_pipeline.sh -- run the full paper pipeline in order, at full scale.
# Not part of the paper code itself; a convenience driver for this run.
set -euo pipefail
cd "$(dirname "$0")"

echo "=== [$(date)] export_empirical_networks.R ==="
(cd networks && Rscript export_empirical_networks.R)

echo "=== [$(date)] 01_ews_sentinel_sweep.jl ==="
julia --project=. -t 15 simulate/01_ews_sentinel_sweep.jl

echo "=== [$(date)] 02_tipmoc_export_csv.jl (mult=1) ==="
julia --project=. simulate/02_tipmoc_export_csv.jl

echo "=== [$(date)] 02_tipmoc_export_csv.jl (mult=2) ==="
julia --project=. simulate/02_tipmoc_export_csv.jl 2

echo "=== [$(date)] 02_tipmoc_detect.py (mult=1) ==="
python3 simulate/02_tipmoc_detect.py

echo "=== [$(date)] 02_tipmoc_detect.py (mult=2) ==="
python3 simulate/02_tipmoc_detect.py 2

echo "=== [$(date)] table1_networks_and_Kc.jl ==="
julia --project=. tables/table1_networks_and_Kc.jl

echo "=== [$(date)] table2_tau_comparison.jl ==="
julia --project=. tables/table2_tau_comparison.jl

echo "=== [$(date)] table3_tipmoc_detection.jl ==="
julia --project=. tables/table3_tipmoc_detection.jl

echo "=== [$(date)] fig_ews_all_kc.jl ==="
julia --project=. figures/fig_ews_all_kc.jl

echo "=== [$(date)] fig_ews_all_kc_appendix.jl ==="
julia --project=. figures/fig_ews_all_kc_appendix.jl

echo "=== [$(date)] fig_ews_sentinel.jl ==="
julia --project=. figures/fig_ews_sentinel.jl

echo "=== [$(date)] fig_ews_sentinel_appendix.jl ==="
julia --project=. figures/fig_ews_sentinel_appendix.jl

echo "=== [$(date)] PIPELINE COMPLETE ==="
