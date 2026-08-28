### Simulation code for "Early warning signals for synchronization transitions from partial observations"

> **Citation:** > Y. Kato and N. Masuda. (2026). *Early warning signals for synchronization transitions from partial observations*. [arXiv:26**.*****](https://arXiv.org/abs/26**.*****)

- The simulation environment (Julia, Python, and R) and necessary packages can be reproduced by Docker. 
- For the brief tutorial of Docker and Apptainer, see our previous GitHub repository: https://github.com/yuu-kato/EMC_kuramoto_2025. 

### Empirical networks
| Directory | File | Description |
| :--- | :--- | :--- |
| `networks/` | `export_empirical_networks.R` | Exports empirical network data from [networkdata](https://schochastics.github.io/networkdata/). |

### Simulation
| File | Description | Results (stored in `cache/`)|
| :--- | :--- | :--- |
| `core/kuramoto_core.jl` | Core functions |
| `simulate/`<br>`01_ews_sentinel_sweep.jl` | Determine critical coupling $K_\mathrm{c}$. <br> Calculate EWSs for $n = \lfloor N \rfloor$. <br> Calculate EWSs for $n = 2\lfloor N \rfloor$. | `kc_summary.tsv` <br> `bundle_ews_Kc20.jls` <br> `bundle_ews_Kc20_mult2.jls` |
|`simulate/`<br>`02_tipmoc_export_csv.jl`|Export results of EWS simulation to .csv files. | `tipmoc_input_Kc20.csv` <br> `tipmoc_input_Kc20_mult2.csv` <br> `tipmoc_meta_Kc20.csv` <br> `tipmoc_meta_Kc20_mult2.csv` |
|`simulate/`<br>`02_tipmoc_detect.py`| Execute TIPMOC detection. | `tipmoc_results_Kc20.csv` <br> `tipmoc_results_Kc20_mult2.csv` |

**Note:** `02_tipmoc_detect.py` uses `delta_aicc.py` and `fit_powerlaw.py` from a local `tipmoc_lib/` directory. These files are identical to those in the [TIPMOC repository](https://github.com/naokimas/tipmoc) from [N. Masuda. (2026)](https://doi.org/10.1093/pnasnexus/pgag126) and are therefore not included in this repository.

### Figures & tables
| File | Produced materials |
| :--- | :--- |
| `figures/fig_ews_all_kc.jl` | EWSs for three representative networks under full observation (S1) (Figure 2).|
| `figures/fig_ews_sentinel.jl` | EWSs for three representative networks under partial observation (S2-S6) (Figure 3).|
| `figures/fig_ews_all_kc_appendix.jl` | EWSs for the other six networks under full observation (S1) (Figure 4).|
| `figures/fig_ews_sentinel_appendix.jl` | EWSs for the other six networks under partial observation (S2-S6) (Figure 5).|
| `tables/table1_networks_and_Kc.jl` | Summary of networks and their critical coupling strengths (Table 1).|
| `tables/table2_tau_comparison.jl` | Kendall's $\tau$ for each combination of EWSs and sentinel sets (Table 2).|
| `tables/table3_tipmoc_detection.jl` | TIPMOC detection rates for each combination of EWSs and sentinel sets (Table 3).|

## Running the pipeline
`run_pipeline.sh` runs every step below in order.

**If you have any questions when running the codes, please contact Yusuke Kato.**