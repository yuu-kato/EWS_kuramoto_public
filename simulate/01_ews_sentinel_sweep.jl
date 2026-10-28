# simulate/01_ews_sentinel_sweep.jl
#
# Critical-coupling determination + main EWS/sentinel-selection driver. 
# Here, `for (i, topo) in enumerate(topos)` builds each
# topology once, and each of its NFREQ=10 frequency realizations builds its
# NetCase exactly once, shared by both the K_c-determination sweep and the
# EWS sweep that follows it in the same job.
#
# For each (network, frequency realization) pair, in one job:
#   stage A : chain an up-sweep K0=0.01 -> Kmax=3.0 (step 0.01; `Kmax_of`) over
#     NRUN_PER_FREQ=10 independent Euler-Maruyama noise realizations,
#     tracking only <R>_t(K) -- no variance of any kind. K_c is the first K
#     at which the run-averaged <R>_t(K) has covered a fraction FRAC of the
#     way from its floor R_0 to full synchronization (R=1)
#     (core/kuramoto_core.jl's `Kc_threshold`). FRAC is the optional first
#     command-line argument (default 0.2, the paper's Kc20); 0.1 and 0.3
#     (Kc10, Kc30) are run as robustness checks of the K_c definition. The
#     resulting K_c is the sole anchor used throughout the rest of the
#     pipeline. The job errors out if any realization's K_c is not reached
#     before Kmax.
#   stage B : anchored to that realization's own K_c, chain a narrower
#     up-sweep K0=0.01 -> 2*K_c for 10 further noise realizations,
#     computing every EWS for all six of the paper's sentinel sets together,
#     from the very same simulated trajectories:
#       S1  all N nodes (full observation)
#       S2  top-n by an EWS-specific criterion (c_i for the order-parameter
#           signals, Var_t delta_theta_i for the microscopic EWS) at K=0.9 K_c
#       S3  same criterion as S2, at K=0.5 K_c
#       S4  top-n by degree (hub nodes)
#       S5  bottom-n by degree (peripheral nodes)
#       S6  n nodes chosen uniformly at random, one independent draw per
#           simulation, the same draw used for all three signals
#     S2 and S3 each use a *signal-dependent* node selection -- ranked by c_i
#     for the local order parameter/collective EWS, or by Var_t delta_theta_i
#     for the microscopic EWS -- i.e. each is really two distinct node sets
#     sharing one K-window. To keep that distinction explicit (rather than
#     silently encoding it as a numeric column index), we track eight
#     underlying node sets per mult, labelled S1, S2-1, S2-2, S3-1, S3-2, S4,
#     S5, S6:
#       S1    all N nodes (full observation)
#       S2-1  top-n by c_i               @ K=0.9Kc  -- R_S and Var_t R_S under S2
#       S2-2  top-n by Var_t delta_theta @ K=0.9Kc  -- microscopic EWS under S2
#       S3-1  top-n by c_i               @ K=0.5Kc  -- R_S and Var_t R_S under S3
#       S3-2  top-n by Var_t delta_theta @ K=0.5Kc  -- microscopic EWS under S3
#       S4    top-n by degree (hub nodes) -- fixed per network, not signal-dependent
#       S5    bottom-n by degree (peripheral nodes) -- fixed per network
#       S6    n nodes chosen uniformly at random -- redrawn every simulation
#     Figures/tables read off the S2-1/S2-2 (resp. S3-1/S3-2) column
#     appropriate to the signal being reported, and S1/S4/S5/S6 directly --
#     see `newset_column` in simulate/02_tipmoc_export_csv.jl, tables/, and
#     figures/ (an identical mapping defined independently in each).
#
# We evaluate two sentinel-count choices, n = mult * floor(ln N) for
# mult in MULTS = [1, 2] (core/kuramoto_core.jl's `n_sentinels`): mult=1 is
# the paper's main definition, mult=2 is the auxiliary "doubled sentinel set
# size" robustness check. Both are computed from the *same* trajectories in
# stage B -- the integration itself (advance_one_K/node_stats_one_K/
# record_one_K) does not depend on sentinel set size at all, only the top-n
# selection and per-set EWS bookkeeping do, which is cheap relative to the
# O(N) or O(E) drift computation performed at every timestep. Re-running the
# whole integration a second time for mult=2 would therefore roughly double
# stage B's cost for no benefit; instead we track 2*8=16 node sets (8 per
# mult) within one pass2 call and split the result into two output bundles
# after the fact.
#
# The hub (S4) and peripheral (S5) sets are a fixed topological property of
# the network, computed once per (network, mult) before the main loop (ties
# in degree are
# broken uniformly at random via a reproducible per-(network,mult,set) seed).
# The complete graph's adjacency is not stored (core/kuramoto_core.jl:
# "unused (all-zero) for :complete"; every node has the same degree N-1 by
# construction), so S4 and S5 there are ill-defined -- both fall back to an
# arbitrary same-size subset, not a meaningful hub/periphery contrast. S6's
# random set is redrawn independently in every (realization, run) *and* per
# mult. Kendall's tau of each EWS vs K is computed over K_0 <= K <= K_c(r), i.e.
# against that realization's own K_c.
#
# Threading flattens over all 9*NFREQ=90 (network, realization) pairs
# (`Threads.@threads`); each job then runs its 10 wide-sweep + 10
# narrow-sweep simulations sequentially.
#
# Run:
#   julia --project=. -t 8 simulate/01_ews_sentinel_sweep.jl        # FRAC = 0.2 
#   julia --project=. -t 8 simulate/01_ews_sentinel_sweep.jl 0.1    # Kc10 robustness check
#   julia --project=. -t 8 simulate/01_ews_sentinel_sweep.jl 0.3    # Kc30 robustness check
#
# Writes under cache/ (names tagged by KCTAG = Kc10/Kc20/Kc30; for the paper's
# FRAC=0.2 the names below are unchanged, other fractions append _Kc10/_Kc30
# to kc_summary and use bundle_ews_Kc10*/bundle_ews_Kc30*):
#   kc_summary.tsv            one row per network: N, E, <k>, mean/sd
#                              of K_c over the 10 realizations (consumed by
#                              tables/table1_networks_and_Kc.jl).
#   bundle_ews_Kc20.jls (mult=1) and bundle_ews_Kc20_mult2.jls (mult=2):
#     `meta` (per-network, realization- and mult-independent: name/label/key/N/E/k),
#     `realization_meta` (per-network-per-realization Kc/Ks/K09/Kmax,
#     identical in both files since it doesn't depend on mult), and
#     `results_allfreq[i][r]` (a length-NRUN_PER_FREQ vector of run records
#     for network i, frequency realization r; `RS`/`RSvar`/`SM` and
#     `tau_RS`/`tau_RSvar`/`tau_micro` are all 8-column, in the
#     S1/S2-1/S2-2/S3-1/S3-2/S4/S5/S6 order above -- see `SET_LABELS`).
# Consumed by figures/fig_ews_sentinel.jl, figures/fig_ews_all_kc.jl,
# simulate/02_tipmoc_export_csv.jl, tables/table3_tipmoc_detection.jl, and
# tables/table2_tau_comparison.jl -- all of which use realization 1 for
# plotted trajectories, but the full 10x10 pool for Kendall-tau summary
# statistics.

include(joinpath(@__DIR__, "..", "core", "kuramoto_core.jl"))
using Serialization, Printf, Statistics, Random

const HERE  = @__DIR__
const CACHE = joinpath(HERE, "..", "cache"); mkpath(CACHE)
const P     = SimParams()
const REGISTRY = network_registry(100)
"Column labels for the eight tracked node sets (per mult) -- see header note
for the S2-1/S2-2 and S3-1/S3-2 signal-dependent-selection explanation."
const SET_LABELS = ["S1", "S2-1", "S2-2", "S3-1", "S3-2", "S4", "S5", "S6"]
const NSETPERMULT = length(SET_LABELS)   # 8
const MULTS = [1, 2]
const FRAC  = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 0.2   # K_c threshold fraction; 0.2 is the paper's
const KCTAG = "Kc" * string(round(Int, 100 * FRAC))               # "Kc10", "Kc20", "Kc30"
const KCSUF = KCTAG == "Kc20" ? "" : "_" * KCTAG                  # keeps kc_summary.tsv unchanged for 0.2

"Sweep range for the K_c-determination grid, same for all networks."
Kmax_of(kind) = 3.0

outfile(mult) = joinpath(CACHE, mult == 1 ? "bundle_ews_$(KCTAG).jls" : "bundle_ews_$(KCTAG)_mult$(mult).jls")

println("Threads = ", Threads.nthreads(), "; anchor = ", KCTAG, "; sentinel mults = ", MULTS,
    "; networks = ", length(REGISTRY), "; freq realizations = ", NFREQ, "; runs/realization = ", NRUN_PER_FREQ)

"Degree vector for `topo`. Complete graph's adjacency is not stored (all
nodes tied at degree N-1 by construction; see module docstring above)."
function degrees_of(topo)
    topo.name == :complete && return fill(Float64(topo.nnodes - 1), topo.nnodes)
    vec(sum(topo.adjacency; dims = 2))
end

"Pass 1: chain K0 -> K09, capturing per-node (vd, ci) at Ks05 and at K09 (all
N nodes -- both mults' top-n selections are truncations of this same output)."
function pass1(net, Ks05, K09, Kgrid, seed)
    rng = MersenneTwister(seed)
    u = 2pi .* rand(rng, net.nnodes)
    vd05 = Float64[]; ci05 = Float64[]; vd09 = Float64[]; ci09 = Float64[]
    for K in Kgrid
        if abs(K - Ks05) < 1e-9
            u, _, vd05, ci05 = node_stats_one_K(K, u, net, P, rng)
        elseif abs(K - K09) < 1e-9
            u, _, vd09, ci09 = node_stats_one_K(K, u, net, P, rng)
            return vd05, ci05, vd09, ci09
        else
            u = advance_one_K(K, u, net, P, rng)
        end
    end
    error("pass1: K09 was never reached on the grid")
end

"Pass 2: full sweep K0 -> 2Kc, tracking every set's three EWS plus the global reference."
function pass2(net, sets, Kgrid, seed)
    rng = MersenneTwister(seed)
    u = 2pi .* rand(rng, net.nnodes)
    nK = length(Kgrid); nsets = length(sets)
    R = zeros(nK); Rvar = zeros(nK); global_micro = zeros(nK)
    RS = zeros(nK, nsets); RSvar = zeros(nK, nsets); set_micro = zeros(nK, nsets)
    for (ik, K) in enumerate(Kgrid)
        u, rm, rv, gm, rsm, rsv, sm = record_one_K(K, u, net, sets, P, rng)
        R[ik] = rm; Rvar[ik] = rv; global_micro[ik] = gm
        @views RS[ik, :] .= rsm; @views RSvar[ik, :] .= rsv; @views set_micro[ik, :] .= sm
    end
    (; R, Rvar, global_micro, RS, RSvar, set_micro)
end

topos = [spec[1]() for spec in REGISTRY]
keys_ = [spec[2] for spec in REGISTRY]
labels = [spec[3] for spec in REGISTRY]
kinds = [spec[4] for spec in REGISTRY]

# ----- degree-based sets S4 (hub) / S5 (periphery): fixed per network, per mult -----
hub_sets = Dict{Tuple{String,Int},Vector{Int}}()
periph_sets = Dict{Tuple{String,Int},Vector{Int}}()
meta = NamedTuple[]
degree_tie_rng = MersenneTwister(20260806_05)
for (i, topo) in enumerate(topos)
    deg = degrees_of(topo)
    for mult in MULTS
        n = n_sentinels(topo; mult = mult)
        hub_sets[(keys_[i], mult)] = select_top(deg, n; rng = degree_tie_rng)
        periph_sets[(keys_[i], mult)] = select_top(-deg, n; rng = degree_tie_rng)
    end
    push!(meta, (; name = topo.name, label = topo.label, key = keys_[i], N = topo.nnodes,
                    E = topo.nedges, k = topo.mean_degree))
    if topo.name == :complete
        println("NOTE: ", topo.label, " -- all degrees tied at N-1; S4/S5 are each an independent random n-node subset (ties broken uniformly at random), not a real hub/periphery contrast.")
    end
end

# ----- main loop: one job per (network, frequency realization); builds its NetCase
# once, determines that realization's own Kc (stage A), then computes every EWS
# anchored to it (stage B) -----
Kc_all = [Vector{Float64}(undef, NFREQ) for _ in topos]
realization_meta = [Vector{NamedTuple}(undef, NFREQ) for _ in topos]
results_allfreq_by_mult = Dict(mult => [Vector{Vector{Any}}(undef, NFREQ) for _ in topos] for mult in MULTS)
for mult in MULTS, i in eachindex(topos), r in 1:NFREQ
    results_allfreq_by_mult[mult][i][r] = Vector{Any}(undef, NRUN_PER_FREQ)
end

jobs = [(i, r) for i in eachindex(topos) for r in 1:NFREQ]
t0 = time()
Threads.@threads for j in 1:length(jobs)
    i, r = jobs[j]
    topo = topos[i]
    net = realize(topo, r)   # built once, shared by stage A and stage B below

    # ---- stage A: wide up-sweep for this realization's own Kc (was 01) ----
    Kgrid_wide = collect(0.01:0.01:Kmax_of(kinds[i]))
    Rmean_wide = Matrix{Float64}(undef, length(Kgrid_wide), NRUN_PER_FREQ)
    for run in 1:NRUN_PER_FREQ
        seed = 20260806_100 + 1000 * run + 37 * r + i
        @views Rmean_wide[:, run] .= sweep_global(net, Kgrid_wide, P; seed = seed)
    end
    Kc_r = Kc_threshold(vec(mean(Rmean_wide; dims = 2)), Kgrid_wide, FRAC)
    Kc_r < Kmax_of(kinds[i]) - 1e-9 || error("Kc never crossed for ", labels[i], " realization ", r, "; raise Kmax_of")
    Kc_all[i][r] = Kc_r

    Kc_snap = snap_to_grid(Kc_r)
    Ks05_r  = snap_to_grid(0.5 * Kc_snap)
    K09_r   = snap_to_grid(0.9 * Kc_snap)
    Kmax_r  = snap_to_grid(2.0 * Kc_snap)
    realization_meta[i][r] = (; Kc = Kc_snap, Ks = Ks05_r, K09 = K09_r, Kmax = Kmax_r)

    # ---- stage B: EWS + sentinel sweep anchored to Kc_snap, same NetCase (was 02) ----
    Kgrid = collect(0.01:0.01:Kmax_r)
    sel_rng = MersenneTwister(20260806_500 + 37 * r + i)
    for run in 1:NRUN_PER_FREQ
        seed_sel = 20260806_600 + 1000 * run + 37 * r + i   # pass1: ranking trajectory for S2/S3 only
        seed     = 20260806_000 + 1000 * run + 37 * r + i   # pass2: evaluation trajectory (unchanged)

        vd05, ci05, vd09, ci09 = pass1(net, Ks05_r, K09_r, Kgrid, seed_sel)

        all_sets = Vector{Vector{Int}}()
        for (mi, mult) in enumerate(MULTS)
            n = n_sentinels(topo; mult = mult)
            s6_set = only(random_sets(net.nnodes, n, 1; rng = sel_rng))
            append!(all_sets, [collect(1:net.nnodes),
                                select_top(ci09, n; rng = sel_rng), select_top(vd09, n; rng = sel_rng),
                                select_top(ci05, n; rng = sel_rng), select_top(vd05, n; rng = sel_rng),
                                hub_sets[(keys_[i], mult)], periph_sets[(keys_[i], mult)], s6_set])   # S1, S2-1, S2-2, S3-1, S3-2, S4, S5, S6
        end

        p2 = pass2(net, all_sets, Kgrid, seed)
        below_Kc = Kgrid .<= Kc_snap + 1e-9
        Km = Kgrid[below_Kc]
        tau_of(M, cols) = [kendall_tau(M[below_Kc, s], Km) for s in cols]

        for (mi, mult) in enumerate(MULTS)
            cols = ((mi - 1) * NSETPERMULT + 1):(mi * NSETPERMULT)
            RS = p2.RS[:, cols]; RSvar = p2.RSvar[:, cols]; SM = p2.set_micro[:, cols]
            tau_RS = tau_of(p2.RS, cols); tau_RSvar = tau_of(p2.RSvar, cols); tau_micro = tau_of(p2.set_micro, cols)
            results_allfreq_by_mult[mult][i][r][run] = (; run, sets = all_sets[cols], Kv = Kgrid,
                R = p2.R, Rvar = p2.Rvar, gm = p2.global_micro,
                RS = RS, RSvar = RSvar, SM = SM,
                tau_RS = tau_RS, tau_RSvar = tau_RSvar, tau_micro = tau_micro)
        end
    end
end
@printf("%d networks x %d freq realizations x %d runs done in %.1f s\n",
    length(topos), NFREQ, NRUN_PER_FREQ, time() - t0)

# ----- Kc summary -----
summary_rows = Vector{NamedTuple}(undef, length(topos))
for (i, topo) in enumerate(topos)
    Kc = Kc_all[i]
    summary_rows[i] = (; key = keys_[i], label = labels[i], kind = kinds[i], N = topo.nnodes, E = topo.nedges,
                         kbar = topo.mean_degree,
                         Kc_mean = mean(Kc), Kc_sd = std(Kc))
    @printf("%-16s N=%4d E=%6d <k>=%7.2f  %s=%.3f(%.3f) range %.2f..%.2f\n",
        labels[i], topo.nnodes, topo.nedges, topo.mean_degree, KCTAG, mean(Kc), std(Kc), minimum(Kc), maximum(Kc))
end

open(joinpath(CACHE, "kc_summary$(KCSUF).tsv"), "w") do io
    println(io, "key\tlabel\tkind\tN\tE\tkbar\tKc_mean\tKc_sd")
    for s in summary_rows
        @printf(io, "%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\n",
            s.key, s.label, s.kind, s.N, s.E, s.kbar, s.Kc_mean, s.Kc_sd)
    end
end
println("wrote kc_summary$(KCSUF).tsv to ", CACHE)

# ----- EWS/sentinel bundles -----
for mult in MULTS
    results_allfreq = results_allfreq_by_mult[mult]
    println("--- mult=", mult, " ---")
    for (i, m) in enumerate(meta)
        all_runs = [results_allfreq[i][r][run] for r in 1:NFREQ for run in 1:NRUN_PER_FREQ]
        tau_micro_mean(s) = mean([rr.tau_micro[s] for rr in all_runs])
        tau_varR_mean(s)  = mean([rr.tau_RSvar[s] for rr in all_runs])
        @printf("%-15s  τ_micro S2=%.2f S3=%.2f S4=%.2f S5=%.2f S6=%.2f | τ_VarR S1=%.2f S2=%.2f S3=%.2f  (pooled over %d runs)\n",
            m.label, tau_micro_mean(3), tau_micro_mean(5), tau_micro_mean(6), tau_micro_mean(7), tau_micro_mean(8),
            tau_varR_mean(1), tau_varR_mean(2), tau_varR_mean(4), length(all_runs))
    end

    serialize(outfile(mult),
        (; kc_mode = KCTAG, params = P, nfreq = NFREQ, nruns = NRUN_PER_FREQ, mult = mult, set_labels = SET_LABELS,
           net_order = [n.name for n in topos], meta = meta, realization_meta = realization_meta,
           results_allfreq = results_allfreq))
    println("saved ", basename(outfile(mult)), " to ", CACHE)
end
