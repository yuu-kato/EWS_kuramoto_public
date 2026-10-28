### Simulation code for "Early warning signals for synchronization transitions from partial observations"

> **Citation:** > Y. Kato and N. Masuda. (2026). *Early warning signals for synchronization transitions from partial observations*. [arXiv:2608.28320](https://arXiv.org/abs/2608.28320)

- The simulation environment (Julia, Python, and R) and necessary packages can be reproduced by Docker. 
- For the brief tutorial of Docker and Apptainer, see our previous GitHub repository: https://github.com/yuu-kato/EMC_kuramoto_2025. 

### Empirical networks
| Directory | File | Description |
| :--- | :--- | :--- |
| `networks/` | `export_empirical_networks.R` | Exports empirical network data from [networkdata](https://schochastics.github.io/networkdata/). |

### Definitions of the critical coupling strength
The critical coupling strength $K_\mathrm{c}$ is the first $K$ at which the normalized order parameter $(R - R_0)/(1 - R_0)$ exceeds a threshold fraction $p$. The paper uses the fraction 0.2 (tag `Kc20`, the default of every script). The fractions 0.1 (`Kc10`) and 0.3 (`Kc30`) are used for robustness checks. All scripts below accept the tag as an argument; without it they read and write the `Kc20` files.

| Tag | Threshold fraction | Suffix of the output files |
| :--- | :--- | :--- |
| `Kc20` | 0.2 (paper) | none |
| `Kc10` | 0.1 | `_Kc10` |
| `Kc30` | 0.3 | `_Kc30` |

### Simulation
| File | Arguments | Description | Results (stored in `cache/`)|
| :--- | :--- | :--- | :--- |
| `core/kuramoto_core.jl` | | Core functions | |
| `simulate/`<br>`01_ews_sentinel_sweep.jl` | `[frac]` (default `0.2`) | Determine critical coupling $K_\mathrm{c}$ with threshold fraction `frac`. <br> Calculate EWSs for $n = \lfloor \ln N \rfloor$. <br> Calculate EWSs for $n = 2\lfloor \ln N \rfloor$. | `kc_summary[_Kc10\|_Kc30].tsv` <br> `bundle_ews_<tag>.jls` <br> `bundle_ews_<tag>_mult2.jls` |
|`simulate/`<br>`02_tipmoc_export_csv.jl`| `[mult] [full\|nofull] [tag]` <br> (defaults `1 nofull Kc20`) | Export results of EWS simulation to .csv files. | `tipmoc_input_<tag>[_mult2][_to2Kc].csv` <br> `tipmoc_meta_<tag>[_mult2].csv` |
|`simulate/`<br>`02_tipmoc_detect.py`| `[mult] [full\|nofull] [tag]` <br> (same as above) | Execute TIPMOC detection.  | `tipmoc_results_<tag>[_mult2][_to2Kc].csv` |

With `nofull`, TIPMOC is applied to $K < K_\mathrm{c}$ only. With `full`, TIPMOC is applied to the whole stored sweep $K_0 \le K \le 2K_\mathrm{c}$, and the output files carry the suffix `_to2Kc`. The `full` option is needed only for `Kc20` with `mult = 1`, because `fig_tipmoc_detection_point.jl` reads `tipmoc_results_Kc20_to2Kc.csv`. For every other combination of tag and `mult`, `nofull` is sufficient.

**Note:** `02_tipmoc_detect.py` uses `delta_aicc.py` and `fit_powerlaw.py` from a local `tipmoc_lib/` directory. These files are identical to those in the [TIPMOC repository](https://github.com/naokimas/tipmoc) from [N. Masuda. (2026)](https://doi.org/10.1093/pnasnexus/pgag126) and are therefore not included in this repository.

### Figures & tables
Outputs are written to `draft/`.

| File | Arguments | Produced materials |
| :--- | :--- | :--- |
| `figures/fig_ews_all_kc.jl` | | EWSs for three representative networks under full observation (S1) (Figure 2).|
| `figures/fig_ews_sentinel.jl` | | EWSs for three representative networks under partial observation (S2-S6) (Figure 3).|
| `figures/fig_ews_all_kc_appendix.jl` | | EWSs for the other six networks under full observation (S1) (Figure 4).|
| `figures/fig_ews_sentinel_appendix.jl` | | EWSs for the other six networks under partial observation (S2-S6) (Figure 5).|
| `figures/fig_tipmoc_detection_point.jl` | | Distribution of the TIPMOC detection point $K_\mathrm{det}/K_\mathrm{c}$ over $K_0 \le K \le 2K_\mathrm{c}$ for each EWS and sentinel set (Figure 6).|
| `tables/table1_networks_and_Kc.jl` | `[tag]` | Summary of networks and their critical coupling strengths (Tables I and IV).|
| `tables/table2_tau_comparison.jl` | `[tag]` | Kendall's $\tau$ for each combination of EWSs and sentinel sets (Tables II and V).|
| `tables/table3_tipmoc_detection.jl` | `[tag]` | TIPMOC detection rates before $K_\mathrm{c}$ for each combination of EWSs and sentinel sets (Tables III and VI).|

### Running the pipeline
`run_pipeline.sh` runs every step in order: the sweep for each of the three $K_\mathrm{c}$ definitions, the TIPMOC export and detection for every (tag, mult) combination, the three tables for each tag, and the figures. Individual steps, e.g.

```bash
julia --project=. -t 8 simulate/01_ews_sentinel_sweep.jl 0.1        # sweep, Kc10
julia --project=. simulate/02_tipmoc_export_csv.jl 1 full            # export Kc20, mult 1, K0..2Kc
python3 simulate/02_tipmoc_detect.py 1 full                          # TIPMOC on that export
julia --project=. simulate/02_tipmoc_export_csv.jl 2 nofull Kc30     # export Kc30, mult 2, K < Kc
python3 simulate/02_tipmoc_detect.py 2 nofull Kc30
julia --project=. tables/table2_tau_comparison.jl Kc10               # Table 2 for Kc10
```

**If you have any questions when running the codes, please contact Yusuke Kato.**
