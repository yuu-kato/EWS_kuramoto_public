#!/usr/bin/env bash
# run_pipeline.sh -- run the full paper pipeline in order, at full scale.
# Not part of the paper code itself; a convenience driver for this run.
#
# K_c definitions: the paper's threshold fraction 0.2 (Kc20, default of
# every script) plus the robustness checks 0.1 (Kc10) and 0.3 (Kc30); the
# sweep is run once per fraction and all downstream scripts take the tag.
# TIPMOC: for Kc20 at n = floor(ln N) the detector is run on the whole
# stored sweep K_0..2K_c (`full`, needed for fig_tipmoc_detection_point.jl);
# every other combination uses the K < K_c window. Table 3 reads the full-
# window results where they exist and the short-window results otherwise
# (only detections before K_c count either way).
set -euo pipefail
cd "$(dirname "$0")"

echo "=== [$(date)] export_empirical_networks.R ==="
(cd networks && Rscript export_empirical_networks.R)

# ---- sweeps: one per K_c threshold fraction ----
echo "=== [$(date)] 01_ews_sentinel_sweep.jl (frac=0.2, Kc20) ==="
julia --project=. -t 15 simulate/01_ews_sentinel_sweep.jl

echo "=== [$(date)] 01_ews_sentinel_sweep.jl (frac=0.1, Kc10) ==="
julia --project=. -t 15 simulate/01_ews_sentinel_sweep.jl 0.1

echo "=== [$(date)] 01_ews_sentinel_sweep.jl (frac=0.3, Kc30) ==="
julia --project=. -t 15 simulate/01_ews_sentinel_sweep.jl 0.3

# ---- TIPMOC: Kc20 ----
echo "=== [$(date)] 02_tipmoc_export_csv.jl (Kc20, mult=1, K0..2Kc) ==="
julia --project=. simulate/02_tipmoc_export_csv.jl 1 full

echo "=== [$(date)] 02_tipmoc_detect.py (Kc20, mult=1, K0..2Kc) ==="
python3 simulate/02_tipmoc_detect.py 1 full

echo "=== [$(date)] 02_tipmoc_export_csv.jl (Kc20, mult=2) ==="
julia --project=. simulate/02_tipmoc_export_csv.jl 2

echo "=== [$(date)] 02_tipmoc_detect.py (Kc20, mult=2) ==="
python3 simulate/02_tipmoc_detect.py 2

# ---- TIPMOC: Kc10 ----
echo "=== [$(date)] 02_tipmoc_export_csv.jl (Kc10, mult=1) ==="
julia --project=. simulate/02_tipmoc_export_csv.jl 1 nofull Kc10

echo "=== [$(date)] 02_tipmoc_detect.py (Kc10, mult=1) ==="
python3 simulate/02_tipmoc_detect.py 1 nofull Kc10

echo "=== [$(date)] 02_tipmoc_export_csv.jl (Kc10, mult=2) ==="
julia --project=. simulate/02_tipmoc_export_csv.jl 2 nofull Kc10

echo "=== [$(date)] 02_tipmoc_detect.py (Kc10, mult=2) ==="
python3 simulate/02_tipmoc_detect.py 2 nofull Kc10

# ---- TIPMOC: Kc30 ----
echo "=== [$(date)] 02_tipmoc_export_csv.jl (Kc30, mult=1) ==="
julia --project=. simulate/02_tipmoc_export_csv.jl 1 nofull Kc30

echo "=== [$(date)] 02_tipmoc_detect.py (Kc30, mult=1) ==="
python3 simulate/02_tipmoc_detect.py 1 nofull Kc30

echo "=== [$(date)] 02_tipmoc_export_csv.jl (Kc30, mult=2) ==="
julia --project=. simulate/02_tipmoc_export_csv.jl 2 nofull Kc30

echo "=== [$(date)] 02_tipmoc_detect.py (Kc30, mult=2) ==="
python3 simulate/02_tipmoc_detect.py 2 nofull Kc30

# ---- tables ----
echo "=== [$(date)] table1_networks_and_Kc.jl (Kc20 / Kc10 / Kc30) ==="
julia --project=. tables/table1_networks_and_Kc.jl
julia --project=. tables/table1_networks_and_Kc.jl Kc10
julia --project=. tables/table1_networks_and_Kc.jl Kc30

echo "=== [$(date)] table2_tau_comparison.jl (Kc20 / Kc10 / Kc30) ==="
julia --project=. tables/table2_tau_comparison.jl
julia --project=. tables/table2_tau_comparison.jl Kc10
julia --project=. tables/table2_tau_comparison.jl Kc30

echo "=== [$(date)] table3_tipmoc_detection.jl (Kc20 / Kc10 / Kc30) ==="
julia --project=. tables/table3_tipmoc_detection.jl
julia --project=. tables/table3_tipmoc_detection.jl Kc10
julia --project=. tables/table3_tipmoc_detection.jl Kc30

# ---- figures ----
echo "=== [$(date)] fig_ews_all_kc.jl ==="
julia --project=. figures/fig_ews_all_kc.jl

echo "=== [$(date)] fig_ews_all_kc_appendix.jl ==="
julia --project=. figures/fig_ews_all_kc_appendix.jl

echo "=== [$(date)] fig_ews_sentinel.jl ==="
julia --project=. figures/fig_ews_sentinel.jl

echo "=== [$(date)] fig_ews_sentinel_appendix.jl ==="
julia --project=. figures/fig_ews_sentinel_appendix.jl

echo "=== [$(date)] fig_tipmoc_detection_point.jl ==="
julia --project=. figures/fig_tipmoc_detection_point.jl

echo "=== [$(date)] PIPELINE COMPLETE ==="