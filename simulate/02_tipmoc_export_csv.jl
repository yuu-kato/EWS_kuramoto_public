# simulate/02_tipmoc_export_csv.jl
#
# Export the EWS trajectories to long-format CSV for TIPMOC (main.tex Sec.
# "Detection of the critical transition by TIPMOC"), from
# cache/bundle_ews_<tag>[_mult<mult>].jls (written by
# simulate/01_ews_sentinel_sweep.jl). Detection itself runs in the original,
# unmodified Python/SciPy implementation (simulate/02_tipmoc_detect.py, using
# tipmoc_lib/delta_aicc.py and fit_powerlaw.py directly) rather than a Julia
# re-implementation -- see that script's header note for why.
#
# Default window: K STRICTLY LESS THAN K_c, not K <= K_c. TIPMOC is an
# early-warning detector, so it is only given the pre-transition portion of
# each trajectory, not the K_c grid point itself. This is a different
# (narrower) window than the one used for Kendall's tau in
# simulate/01_ews_sentinel_sweep.jl, which does include K = K_c -- that's a
# summary statistic over "the approach to K_c", not an early-warning
# detector, so including the endpoint there is appropriate.
#
# With the second argument `full`, the whole stored range K_0..2K_c is
# exported instead (the bundle's sweep already ends at 2K_c, so no new
# simulation is needed). This feeds the late-detection analysis: where does
# TIPMOC first fire when it is allowed to run past the transition? Output
# names then carry the suffix _to2Kc so the default files are not touched.
#
# Exported scopes are S1-S6, the six sentinel sets tracked in the bundle (not
# the 8 underlying node-set columns -- see simulate/01_ews_sentinel_sweep.jl's
# header note): for
# signal "R_S"/"VarR" (local order parameter / collective EWS), S2 and S3
# resolve to the bundle's S2-1/S3-1 (c_i-ranked) columns; for signal "micro"
# (microscopic EWS), they resolve to S2-2/S3-2 (Var_t delta_theta-ranked).
# S1, S4 (hub), S5 (periphery), and S6 (random) resolve to a single column
# regardless of signal -- see `newset_column` below, identical mapping to
# figures/fig_ews_sentinel.jl and tables/table2_tau_comparison.jl.
#
# All NFREQ=10 frequency realizations x NRUN_PER_FREQ=10 runs per network
# (all 100 trajectories) are exported. `run` is a *global* run index 1..100
# (`(freq-1)*NRUN_PER_FREQ + run_within_freq`), since different frequency
# realizations have different anchor K_c values and therefore need their own
# entries in tipmoc_meta -- see that file's schema below.
#
# Optional positional command-line arguments, in this order:
#   1. `mult`  sentinel-count multiplier (default 1); must match the bundle.
#   2. `full`  export K_0..2K_c instead of K < K_c (default: K < K_c; pass any
#              other word, e.g. `nofull`, as a placeholder if argument 3 is
#              needed without it).
#   3. tag     K_c-definition family of the bundle to export: Kc20 (default,
#              the paper's threshold fraction 0.2), Kc10 or Kc30.
#
# Run (after simulate/01_ews_sentinel_sweep.jl with the matching fraction):
#   julia --project=. simulate/02_tipmoc_export_csv.jl                 # mult 1, K<Kc, Kc20
#   julia --project=. simulate/02_tipmoc_export_csv.jl 2               # mult 2
#   julia --project=. simulate/02_tipmoc_export_csv.jl 1 full          # K0..2Kc
#   julia --project=. simulate/02_tipmoc_export_csv.jl 1 full Kc10     # K0..2Kc, Kc10 bundle
#
# Writes cache/tipmoc_input_<tag>[_mult<mult>][_to2Kc].csv (columns: net,scope,
# signal,run,k,value; run is the global 1..100 index) and
# cache/tipmoc_meta_<tag>[_mult<mult>].csv (columns: net,run,Kc -- one row per
# (net, global run index), since Kc depends on the frequency realization; the
# same meta file serves both windows), consumed by simulate/02_tipmoc_detect.py.

include(joinpath(@__DIR__, "..", "core", "kuramoto_core.jl"))
using Serialization: deserialize
using Printf

const HERE  = @__DIR__
const CACHE = joinpath(HERE, "..", "cache")
const MULT = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1
const TAG = MULT == 1 ? "" : "_mult$(MULT)"
const FULL = length(ARGS) >= 2 && ARGS[2] == "full"   # export K0..2Kc instead of K < Kc
const RTAG = FULL ? "_to2Kc" : ""
const KCTAG = length(ARGS) >= 3 ? ARGS[3] : "Kc20"   # which K_c-definition bundle to export: Kc10, Kc20, Kc30
const B = deserialize(joinpath(CACHE, "bundle_ews_$(KCTAG)$(TAG).jls"))
const SCOPES = ["S1", "S2", "S3", "S4", "S5", "S6"]
const SIGNALS = [("R_S", :RS), ("VarR", :RSvar), ("micro", :SM)]

"Bundle column (1..8, S1/S2-1/S2-2/S3-1/S3-2/S4/S5/S6 -- see simulate/02's
header note) for sentinel set `scope` (S1..S6) and `signal` ('R_S', 'VarR',
or 'micro')."
function newset_column(scope, signal)
    scope == "S1" && return 1
    scope == "S2" && return signal == "micro" ? 3 : 2
    scope == "S3" && return signal == "micro" ? 5 : 4
    scope == "S4" && return 6   # hub
    scope == "S5" && return 7   # periphery
    scope == "S6" && return 8   # random
    error("scope must be one of $SCOPES, got $scope")
end

global_run(freq, run) = (freq - 1) * B.nruns + run   # 1..NFREQ*NRUN_PER_FREQ, e.g. 1..100

open(joinpath(CACHE, "tipmoc_input_$(KCTAG)$(TAG)$(RTAG).csv"), "w") do io
    println(io, "net,scope,signal,run,k,value")
    for (i, m) in enumerate(B.meta)
        for freq in 1:B.nfreq
            Kc = B.realization_meta[i][freq].Kc
            for rr in B.results_allfreq[i][freq]
                grun = global_run(freq, rr.run)
                for scope in SCOPES
                    for (signame, fld) in SIGNALS
                        col = newset_column(scope, signame)
                        values = getfield(rr, fld)
                        for ik in eachindex(rr.Kv)
                            FULL || rr.Kv[ik] < Kc - 1e-9 || continue   # strictly K < K_c unless FULL (then K0..2Kc)
                            @printf(io, "%s,%s,%s,%d,%.4f,%.10g\n", m.name, scope, signame, grun, rr.Kv[ik], values[ik, col])
                        end
                    end
                end
            end
        end
    end
end
open(joinpath(CACHE, "tipmoc_meta_$(KCTAG)$(TAG).csv"), "w") do io
    println(io, "net,run,Kc")
    for (i, m) in enumerate(B.meta)
        for freq in 1:B.nfreq
            Kc = B.realization_meta[i][freq].Kc
            for run in 1:B.nruns
                @printf(io, "%s,%d,%.4f\n", m.name, global_run(freq, run), Kc)
            end
        end
    end
end
println("wrote tipmoc_input_$(KCTAG)$(TAG)$(RTAG).csv and tipmoc_meta_$(KCTAG)$(TAG).csv to ", CACHE,
    " (", length(B.meta), " networks x ", B.nfreq * B.nruns, " trajectories each, ", length(SCOPES), " scopes x ", length(SIGNALS), " signals)")
