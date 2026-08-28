"""simulate/02_tipmoc_detect.py

Apply TIPMOC (https://github.com/naokimas/tipmoc; vendored, unmodified, at
../tipmoc_lib -- see main.tex Sec. "Detection of the critical transition
by TIPMOC") to every (network, sentinel scope, EWS signal, run) trajectory in
cache/tipmoc_input_Kc20[_mult<mult>].csv (written by
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

The (optional) command-line argument is the sentinel-count multiplier `mult`
(default 1) and must match whichever tipmoc_input_Kc20[_mult<mult>].csv you
want to process (see simulate/02_tipmoc_export_csv.jl).

Run (after simulate/02_tipmoc_export_csv.jl):
    python3 simulate/02_tipmoc_detect.py
    python3 simulate/02_tipmoc_detect.py 2

Writes cache/tipmoc_results_Kc20[_mult<mult>].csv (columns: net,scope,
signal,run,n,Kc,detected,u_det,hatuc_det,gamma_det,corr_det) -- one row per
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
        return {**base, "detected": False, "u_det": np.nan, "hatuc_det": np.nan, "gamma_det": np.nan, "corr_det": np.nan}

    n_windows = n - I0 + 1
    dA = np.full(n_windows, np.nan)
    hatuc = dA.copy(); gamma = dA.copy(); corr = dA.copy(); u_end = dA.copy()
    for k, i in enumerate(range(I0, n + 1)):
        d, h, g_, c = _delta(u[:i], v[:i])
        dA[k], hatuc[k], gamma[k], corr[k], u_end[k] = d, h, g_, c, u[i - 1]

    detected = False
    u_det = hatuc_det = gamma_det = corr_det = np.nan
    for j in range(len(dA) - 2):
        window = dA[j:j + 3]
        if np.all(np.isfinite(window)) and np.all(window <= THRESH):
            detected = True
            u_det, hatuc_det, gamma_det, corr_det = u_end[j + 2], hatuc[j + 2], gamma[j + 2], corr[j + 2]
            break
    return {**base, "detected": detected, "u_det": u_det, "hatuc_det": hatuc_det, "gamma_det": gamma_det, "corr_det": corr_det}


def main():
    cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cache")
    df = pd.read_csv(os.path.join(cache, f"tipmoc_input_Kc20{TAG}.csv"))
    meta = pd.read_csv(os.path.join(cache, f"tipmoc_meta_Kc20{TAG}.csv")).set_index(["net", "run"])["Kc"].to_dict()
    groups = [((net, scope, signal, run), g, meta[(net, run)])
              for (net, scope, signal, run), g in df.groupby(["net", "scope", "signal", "run"])]
    print(len(groups), "groups, mult =", MULT, flush=True)
    # "spawn" avoids inheriting a forked, possibly-inconsistent BLAS thread pool.
    # Leave a couple of cores free rather than saturating the machine.
    nproc = max(1, (os.cpu_count() or 1) - 2)
    ctx = mp.get_context("spawn")
    rows = []
    with ctx.Pool(processes=nproc) as pool:
        for i, row in enumerate(pool.imap_unordered(run_group, groups, chunksize=4), 1):
            rows.append(row)
            if i % 50 == 0 or i == len(groups):
                print(f"  {i}/{len(groups)} groups done", flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(cache, f"tipmoc_results_Kc20{TAG}.csv"), index=False)
    print(f"wrote tipmoc_results_Kc20{TAG}.csv", out.shape, "detection rate", round(out["detected"].mean(), 3))


if __name__ == "__main__":
    main()
