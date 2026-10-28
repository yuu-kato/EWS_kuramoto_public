"""simulate/02_tipmoc_detect.py

Apply TIPMOC (https://github.com/naokimas/tipmoc; vendored, unmodified, at
../tipmoc_lib -- see main.tex Sec. "Detection of the critical transition
by TIPMOC") to every (network, sentinel scope, EWS signal, run) trajectory in
cache/tipmoc_input_<tag>[_mult<mult>][_to2Kc].csv (written by
simulate/02_tipmoc_export_csv.jl). This calls the original, unmodified
delta_AICc / fit_powerlaw (tipmoc_lib/delta_aicc.py, fit_powerlaw.py --
scipy.optimize.differential_evolution with its real defaults: strategy=
best1bin, dithered mutation factor, and a final polish=True local L-BFGS-B
refinement) directly, rather than a Julia re-implementation of them. 

TIPMOC scans an expanding window of the trajectory and, at each window,
compares a power-law divergence fit against a linear null model via the
corrected Akaike information criterion (AICc). A transition is declared when
three consecutive windows all have Delta-AICc <= THRESH (decisively favoring
the power-law/tipping model); I0 is the minimum window length before the
first evaluation. These parameters (I0=8, THRESH=-10, 3 consecutive windows)
match TIPMOC's own reference usage and must not be changed independently of
the upstream method.

The scan stops at the first window that completes such a triple (early
stopping). This is the same decision rule as evaluating every window and
then searching for the first triple, but avoids fitting windows that can no
longer change the outcome; u_det is the K at the right edge of that window,
i.e. the first K at which TIPMOC declares the transition. With the K_0..2K_c
input (see below) u_det may therefore lie at or beyond K_c; tables/
table3_tipmoc_detection.jl counts only u_det < K_c as a detection, so the
reported rate is the same "detected before K_c" statistic for both inputs.

Optional positional command-line arguments, in this order (must match the
export that produced the input):
    1. mult   sentinel-count multiplier (default 1)
    2. full   process the K_0..2K_c export (default: the K < K_c export; pass
              any other word, e.g. nofull, as a placeholder if argument 3 is
              needed without it)
    3. tag    K_c-definition family: Kc20 (default), Kc10 or Kc30

Run (after simulate/02_tipmoc_export_csv.jl with the same arguments):
    python3 simulate/02_tipmoc_detect.py
    python3 simulate/02_tipmoc_detect.py 2
    python3 simulate/02_tipmoc_detect.py 1 full
    python3 simulate/02_tipmoc_detect.py 1 full Kc10

Writes cache/tipmoc_results_<tag>[_mult<mult>][_to2Kc].csv (columns: net,scope,
signal,run,n,Kc,detected,u_det,hatuc_det,gamma_det,corr_det,n_eval; n is the
number of exported K points, n_eval the window length at which the scan
stopped, Kc the realization's anchor copied from the meta file) -- one row per
(network, scope, signal, run) trajectory -- consumed by
tables/table3_tipmoc_detection.jl, which pools over runs into the per-cell
detection rate (%) the table reports.
"""
import os
# Keep each worker process single-threaded: with one worker process per core
# each potentially spawning its own BLAS thread pool, oversubscription slows
# things down far more than it helps.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import numpy as np
import pandas as pd
import multiprocessing as mp

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tipmoc_lib"))
import delta_aicc as dAICc

I0 = 8
THRESH = -10.0
MULT = int(sys.argv[1]) if len(sys.argv) > 1 else 1
TAG = "" if MULT == 1 else f"_mult{MULT}"
FULL = len(sys.argv) > 2 and sys.argv[2] == "full"
RTAG = "_to2Kc" if FULL else ""
KCTAG = sys.argv[3] if len(sys.argv) > 3 else "Kc20"   # Kc10, Kc20, Kc30


def _delta(u, v):
    """Call TIPMOC's delta_AICc and normalize failures to a NaN tuple."""
    try:
        out = dAICc.delta_AICc(u, v)
    except Exception:
        return (np.nan, np.nan, np.nan, np.nan)
    if out is None or out[0] is None or len(out) != 7:
        return (np.nan, np.nan, np.nan, np.nan)
    d, _popt_lin, _a, _b, hatuc, gamma, corr = out
    return (float(d), float(hatuc), float(gamma), float(corr))


def run_group(args):
    (net, scope, signal, run), g, Kc = args
    g = g.sort_values("k")
    u = g["k"].to_numpy(float)
    v = g["value"].to_numpy(float)
    n = len(u)
    base = dict(net=net, scope=scope, signal=signal, run=run, n=n, Kc=Kc)
    if n < I0:
        return {**base, "detected": False, "u_det": np.nan, "hatuc_det": np.nan, "gamma_det": np.nan, "corr_det": np.nan, "n_eval": 0}
    
    detected = False
    u_det = hatuc_det = gamma_det = corr_det = np.nan
    last3 = []                                    # Delta-AICc of the three most recent windows
    n_eval = 0
    for i in range(I0, n + 1):
        d, h, g_, c = _delta(u[:i], v[:i])
        n_eval = i
        last3 = (last3 + [d])[-3:]
        if len(last3) == 3 and np.all(np.isfinite(last3)) and np.all(np.array(last3) <= THRESH):
            detected = True
            u_det, hatuc_det, gamma_det, corr_det = u[i - 1], h, g_, c
            break
    return {**base, "detected": detected, "u_det": u_det, "hatuc_det": hatuc_det,
            "gamma_det": gamma_det, "corr_det": corr_det, "n_eval": n_eval}

def main():
    cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cache")
    df = pd.read_csv(os.path.join(cache, f"tipmoc_input_{KCTAG}{TAG}{RTAG}.csv"))
    meta = pd.read_csv(os.path.join(cache, f"tipmoc_meta_{KCTAG}{TAG}.csv")).set_index(["net", "run"])["Kc"].to_dict()
    groups = [((net, scope, signal, run), g, meta[(net, run)])
              for (net, scope, signal, run), g in df.groupby(["net", "scope", "signal", "run"])]
    print(len(groups), "groups, mult =", MULT, "full =", FULL, flush=True)
    # "spawn" avoids inheriting a forked, possibly-inconsistent BLAS thread pool.
    # Leave a couple of cores free rather than saturating the machine.
    nproc = max(1, (os.cpu_count() or 1) - 2)
    ctx = mp.get_context("spawn")
    rows = []
    with ctx.Pool(processes=nproc) as pool:
        for i, row in enumerate(pool.imap_unordered(run_group, groups, chunksize=1), 1):
            rows.append(row)
            if i % 50 == 0 or i == len(groups):
                print(f"  {i}/{len(groups)} groups done", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(cache, f"tipmoc_results_{KCTAG}{TAG}{RTAG}.csv"), index=False)
    print(f"wrote tipmoc_results_{KCTAG}{TAG}{RTAG}.csv", out.shape, "detection rate", round(out["detected"].mean(), 3))


if __name__ == "__main__":
    main()