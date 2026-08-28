# kuramoto_core.jl
#
# Shared simulation and analysis library for "Detecting early warning signals
# for synchronization transitions from partial observations."
#
# This is the reorganized "code for paper" version of core/kuramoto_core.jl.
# Relative to the original it has two substantive changes (both requested for
# the paper-release cleanup):
#
#   1. All temporal-variance estimators are now UNBIASED sample variances
#      rather than the population ("divide by M") estimates the original
#      pipeline used:
#        - Var_t R(t) and Var_t R_S(t) (order-parameter variance, global or
#          over a sentinel set) divide the sum of squared deviations by
#          (M - 1), M = number of sampled time points -- the standard
#          Bessel correction for a mean estimated from the same M samples.
#        - Var_t[delta_theta_i(t)] (the microscopic EWS) divides by (M - 2),
#          not (M - 1), because delta_theta_i's "mean" here is really a
#          two-parameter (slope, intercept) OLS line fit to the M samples,
#          which removes 2 degrees of freedom rather than 1; see `resid_var`
#          below.
#      See `var_unbiased` and `resid_var`.
#
#   2. Sentinel-set numbering: the original pipeline had S1-S4 (full
#      observation, two EWS-ranked sets, and a uniformly-random set) as its
#      main definition, with hub/periphery sets bolted on afterwards as S5/S6
#      in a second, separate sweep script. Here the six main-text sentinel
#      sets are renumbered so the two purely-topological (degree-based) sets
#      sit next to each other and the random set moves last:
#        S1  full observation (all N nodes)
#        S2  top-n by EWS-specific criterion (c_i or Var_t delta_theta_i) at K=0.9 K_c
#        S3  top-n by EWS-specific criterion at K=0.5 K_c
#        S4  top-n by degree (hub nodes)            
#        S5  bottom-n by degree (peripheral nodes)   
#        S6  n nodes chosen uniformly at random      
#      All six are now computed together, from the same simulated
#      trajectories, in simulate/01_ews_sentinel_sweep.jl -- see that file's
#      header note.
#
# Model (main.tex, Eq. 1): noisy Kuramoto oscillators on a network,
#   dtheta_i = [omega_i + (K/<k>) * sum_j A_ij sin(theta_j - theta_i)] dt + sigma dW_i,
# integrated by Euler-Maruyama with dt=0.01, sigma=0.15, for T=300 time units
# (150 burn-in), sampling every 0.1 time units.
#
# Each network's topology (adjacency) is built once and fixed throughout the
# study. Natural frequencies omega_i are deterministic Gaussian quantiles of
# N(0,1) (Sec. "Numerical simulations" of main.tex); UNLIKE topology, the
# node-to-frequency assignment is *not* held fixed -- for every network we
# draw NFREQ=10 independent random assignments (`realize(topology, r)` for
# r in 1:NFREQ), and for each of those we run NRUN_PER_FREQ=10 independent
# noise realizations, for NFREQ*NRUN_PER_FREQ=100 simulations per network in
# total. The empirical critical coupling K_c is therefore also
# assignment-specific: each of the 10 frequency realizations gets its own
# K_c (from its own 10 noise runs), and Kendall's tau / TIPMOC detection are
# evaluated against that realization's own K_c. Figures that show raw
# EWS-vs-K trajectories (rather than a derived summary statistic) use only
# the first frequency realization (r=1) and its second noise run (run=2),
# shown raw and unaveraged, to keep the plots readable -- see
# simulate/01_ews_sentinel_sweep.jl and the figures/ scripts.
#
# Nine networks: three model networks (all-to-all, Erdos-Renyi, scale-free
# gamma=3 via the configuration model), all N=100, and six empirical networks
# (Bat, Highschool boy, House finch, Surfer, Hall, Jazz) loaded from
# networks/nets/<name>_edges.csv and reduced to their largest connected
# component (Sec. "Networks" / Appendix A of main.tex).
#
# This file only defines data structures and functions; it does not run
# anything. It is `include`-d by every script under simulate/, tables/, and
# figures/.

using Statistics, Distributions, Random, SparseArrays, LinearAlgebra
import Graphs   # only used for make_scalefree's configuration-model wiring; kept qualified to avoid clashing with this file's own is_connected

# =============================================================================
# Data structures
# =============================================================================

Base.@kwdef struct SimParams
    noise_sigma::Float64 = 0.15
    dt::Float64          = 0.01
    T::Float64           = 300.0
    burnin::Float64      = 150.0
    sample_dt::Float64   = 0.1
end

"A network's fixed topology (built once), with no frequency assignment yet.
`freq_seed_base` seeds the family of NFREQ independent frequency realizations
drawn from this topology via `realize`."
Base.@kwdef struct NetTopology
    name::Symbol
    label::String
    adjacency::SparseMatrixCSC{Float64,Int}   # unused (all-zero) for :complete; see drift!
    mean_degree::Float64
    seed_used::Int
    nnodes::Int
    nedges::Int
    freq_seed_base::Int
end

"A topology combined with one specific frequency realization; this is what
the simulation routines (drift!, advance_one_K, ...) operate on."
Base.@kwdef struct NetCase
    name::Symbol
    label::String
    adjacency::SparseMatrixCSC{Float64,Int}   # unused (all-zero) for :complete; see drift!
    mean_degree::Float64
    omega::Vector{Float64}
    seed_used::Int
    nnodes::Int
    nedges::Int
end

"Number of independent frequency-to-node assignments drawn per network."
const NFREQ = 10
"Number of independent noise realizations (runs) per frequency assignment."
const NRUN_PER_FREQ = 10

"""
    realize(topo::NetTopology, r::Int) -> NetCase

Attach the `r`-th (1..NFREQ) independent frequency realization to a fixed
topology. Different `r` give different node-to-frequency assignments (via
distinct seeds derived from `topo.freq_seed_base`); the topology itself is
unchanged.
"""
realize(topo::NetTopology, r::Int) =
    NetCase(name = topo.name, label = topo.label, adjacency = topo.adjacency,
            mean_degree = topo.mean_degree, omega = freqs_random(topo.nnodes; seed = topo.freq_seed_base + r),
            seed_used = topo.seed_used, nnodes = topo.nnodes, nedges = topo.nedges)

# =============================================================================
# Natural frequencies
# =============================================================================

"""
    gauss_quantiles(N; m=0.0, s=1.0)

The N deterministic quantiles of a Gaussian(mean=m, sd=s), i.e. omega_j such
that Phi((omega_j - m)/s) = (j - 0.5)/N. Using quantiles rather than random
draws avoids sampling fluctuations in the empirical frequency distribution
g(omega) (Hong et al. 2015; main.tex, Sec. "Numerical simulations").
"""
gauss_quantiles(N; m=0.0, s=1.0) = m .+ s .* quantile.(Normal(), ((1:N) .- 0.5) ./ N)

"""
    freqs_random(N; m=0.0, s=1.0, seed)

Deterministic Gaussian quantiles N(m, s^2), randomly assigned to the N nodes
with a fixed seed. The manuscript sets g(omega) = N(0, 1).
"""
freqs_random(N; m=0.0, s=1.0, seed=20260619) = shuffle(MersenneTwister(seed), gauss_quantiles(N; m=m, s=s))

"""
    Kc_threshold(R::AbstractVector, K::AbstractVector, frac; R0_npts=1)

First K at which R has covered a fraction `frac` of the way from its floor
R0 to full synchronization (R=1): the first K with (R(K)-R0)/(1-R0) > frac.
R0 is read off the first `R0_npts` grid point(s); by default (R0_npts=1)
this is simply R at the first point of the simulation grid, K_0=0.01 (main.tex
Eq. def_kcp). Computed entirely from the <R>_t(K) curve itself -- no
reference to Var_t R(t) or any other externally-derived quantity. Used with
frac=0.2 as `Kc20` (simulate/01_ews_sentinel_sweep.jl), the sole critical-
coupling anchor used throughout this codebase.
"""
function Kc_threshold(R::AbstractVector, K::AbstractVector, frac::Real; R0_npts::Int = 1)
    R0 = mean(@view R[1:R0_npts])
    target = R0 + frac * (1.0 - R0)
    i = findfirst(>(target), R)
    i === nothing ? last(K) : K[i]
end

# =============================================================================
# Graph utilities (connectivity / largest connected component)
# =============================================================================

function is_connected(adj::SparseMatrixCSC)
    N = size(adj, 1)
    N == 0 && return false
    visited = falses(N)
    queue = Int[1]
    visited[1] = true
    rows = rowvals(adj)
    head = 1
    while head <= length(queue)
        j = queue[head]; head += 1
        for r in nzrange(adj, j)
            i = rows[r]
            if !visited[i]
                visited[i] = true
                push!(queue, i)
            end
        end
    end
    all(visited)
end

function sparse_from_edges(edges::Vector{Tuple{Int,Int}}, N)
    rows = Int[]; cols = Int[]
    for (i, j) in edges
        push!(rows, i); push!(cols, j)
        push!(rows, j); push!(cols, i)
    end
    A = sparse(rows, cols, ones(Float64, length(rows)), N, N)
    fill!(nonzeros(A), 1.0)   # collapse any duplicate edges to weight 1
    A
end

nedges_of(adj::SparseMatrixCSC) = Int(nnz(adj) ÷ 2)

"Spacing between successive retry seeds in `make_er`/`make_scalefree`'s
regenerate-until-valid loops, large enough that the two functions' retry
windows (`seed:seed+max_attempts-1`, with default `seed`s only 2 apart) can
never land on the same absolute seed -- while leaving each function's `a=0`
seed (the one actually used in practice; see their docstrings) unchanged,
since a stride only affects `a >= 1`."
const GRAPH_RETRY_STRIDE = 1_000_000

"Node indices of the largest connected component of `adj`."
function largest_cc(adj::SparseMatrixCSC)
    N = size(adj, 1)
    rows = rowvals(adj)
    comp = zeros(Int, N)
    ncomp = 0
    for s in 1:N
        comp[s] == 0 || continue
        ncomp += 1
        queue = Int[s]; comp[s] = ncomp; head = 1
        while head <= length(queue)
            j = queue[head]; head += 1
            for r in nzrange(adj, j)
                i = rows[r]
                if comp[i] == 0
                    comp[i] = ncomp
                    push!(queue, i)
                end
            end
        end
    end
    sizes = zeros(Int, ncomp)
    for c in comp; sizes[c] += 1; end
    findall(==(argmax(sizes)), comp)
end

# =============================================================================
# Network constructors
# =============================================================================

"Complete graph (all-to-all coupling). The adjacency matrix is not stored
(see drift!); mean_degree = N-1 so that the coupling term below reduces to
the standard mean-field Kuramoto form K/N * sum_j sin(theta_j - theta_i)."
make_complete(N; freq_seed_base = 20260619) =
    NetTopology(name = :complete, label = "complete", adjacency = spzeros(Float64, N, N),
                mean_degree = Float64(N - 1), seed_used = 0, nnodes = N,
                nedges = N * (N - 1) ÷ 2, freq_seed_base = freq_seed_base)

"""
Erdos-Renyi G(N,p) graph, regenerated with a new seed until connected. `p`
defaults to `10/(N-1)` so the mean degree is about 10 regardless of N.
"""
function make_er(N; p = 10 / (N - 1), seed = 20260701, freq_seed_base = 20260619, max_attempts = 2000)
    for a in 0:max_attempts-1
        rng = MersenneTwister(seed + GRAPH_RETRY_STRIDE * a)
        edges = Tuple{Int,Int}[]
        for j in 2:N, i in 1:j-1
            rand(rng) < p && push!(edges, (i, j))
        end
        adj = sparse_from_edges(edges, N)
        if is_connected(adj)
            return NetTopology(name = :erdos_renyi, label = "ER", adjacency = adj,
                                mean_degree = mean(vec(sum(adj; dims = 2))), seed_used = seed + GRAPH_RETRY_STRIDE * a,
                                nnodes = N, nedges = nedges_of(adj), freq_seed_base = freq_seed_base)
        end
    end
    error("make_er: failed to find a connected realization within $max_attempts attempts")
end

"""
    make_scalefree(N; gamma=3.0, kmin=4, ...)

Scale-free network with degree exponent `gamma`: a degree sequence is drawn
from the truncated discrete power law P(k) ∝ k^-gamma (kmin <= k <= N-1) and
wired into a simple graph with the *exact* target degree sequence using
Graphs.jl's configuration-model implementation, `Graphs.random_configuration_model`.
Unlike a naive single-pass stub-matching-then-discard construction (which
under-realizes high-degree nodes because rejected self-loop/duplicate stubs
are simply dropped), Graphs.jl recycles rejected stubs into further rounds of
matching until every stub is consumed into a valid edge, so the realized
degree sequence matches the drawn one exactly (up to the `check_graphical`
sequence being realizable as a simple graph at all -- see below).
Regenerated with a new seed until the drawn degree sequence is graphical and
the realized graph is connected.
"""
function make_scalefree(N; gamma = 3.0, kmin = 4, seed = 20260703, freq_seed_base = 20260619, max_attempts = 4000)
    ks = collect(kmin:N-1)
    w = Float64.(ks) .^ (-gamma)
    w ./= sum(w)
    cw = cumsum(w)
    for a in 0:max_attempts-1
        rng = MersenneTwister(seed + GRAPH_RETRY_STRIDE * a)
        draw() = ks[searchsortedfirst(cw, rand(rng))]
        degs = [draw() for _ in 1:N]
        isodd(sum(degs)) && (degs[argmin(degs)] += 1)   # need an even stub count

        g = try
            # check_graphical=true raises ArgumentError (caught below) instead of
            # hanging, since a non-graphical sequence can never be fully matched.
            Graphs.random_configuration_model(N, degs; rng = rng, check_graphical = true)
        catch err
            err isa ArgumentError || rethrow()
            continue   # this draw wasn't graphical; try a fresh degree sequence
        end
        Graphs.is_connected(g) || continue

        edge_list = [(Int(Graphs.src(e)), Int(Graphs.dst(e))) for e in Graphs.edges(g)]
        adj = sparse_from_edges(edge_list, N)
        return NetTopology(name = :scale_free, label = "scale-free", adjacency = adj,
                            mean_degree = mean(vec(sum(adj; dims = 2))), seed_used = seed + GRAPH_RETRY_STRIDE * a,
                            nnodes = N, nedges = nedges_of(adj), freq_seed_base = freq_seed_base)
    end
    error("make_scalefree: failed to find a connected realization within $max_attempts attempts")
end

"""
    load_empirical(name; label, dir, freq_seed)

Load an undirected, unweighted empirical network from `dir/<name>_edges.csv`
(as exported by networks/export_empirical_networks.R) and reduce it to its
largest connected component.
"""
function load_empirical(name::AbstractString; label, dir, freq_seed_base)
    path = joinpath(dir, string(name, "_edges.csv"))
    node_labels = String[]
    node_index = Dict{String,Int}()
    pairs = Set{Tuple{Int,Int}}()
    open(path) do f
        first_line = true
        for ln in eachline(f)
            if first_line
                first_line = false
                continue
            end
            isempty(strip(ln)) && continue
            cols = split(ln, ',')
            a = strip(cols[1], ['"', ' ', '\r', '\n'])
            b = strip(cols[2], ['"', ' ', '\r', '\n'])
            ia = get!(node_index, a) do
                push!(node_labels, a); length(node_labels)
            end
            ib = get!(node_index, b) do
                push!(node_labels, b); length(node_labels)
            end
            ia == ib && continue
            push!(pairs, ia < ib ? (ia, ib) : (ib, ia))
        end
    end
    N = length(node_labels)
    adj = sparse_from_edges(collect(pairs), N)
    keep = largest_cc(adj)
    adj = adj[keep, keep]
    N = length(keep)
    NetTopology(name = Symbol(name), label = label, adjacency = adj,
                mean_degree = mean(vec(sum(adj; dims = 2))), seed_used = 0, nnodes = N,
                nedges = nedges_of(adj), freq_seed_base = freq_seed_base)
end

# =============================================================================
# Network registry (single source of truth for the 9 networks in the paper)
# =============================================================================

const NETS_DIR = joinpath(@__DIR__, "..", "networks", "nets")

"""
    network_registry(N=100)

Returns a Vector of (builder, key, label, kind) tuples for the nine networks
used throughout the paper, in the order of Table I. Each builder returns a
`NetTopology` (built once; combine with `realize(topo, r)` to get the r-th of
NFREQ frequency realizations as a `NetCase`). `kind` is :model or :empirical
(only used to pick the K sweep range in simulate/01; see Kmax_of therein).
`N` only applies to the three model networks (all-to-all, ER, scale-free);
empirical networks use their own node count. Each of the nine networks is
given its own `freq_seed_base`, spaced 1000 apart (rather than 1 apart) so
that no two networks' `realize(topo, r)` seed windows (`freq_seed_base+1` ..
`freq_seed_base+NFREQ`) can overlap -- with adjacent bases the two windows
would share NFREQ-1 seeds, which for same-N networks (the three model
networks all use N=100) means literally identical `omega` assignments (same
seed + same `gauss_quantiles(N)` input to `shuffle`).
"""
network_registry(N = 100) = [
    (() -> make_complete(N; freq_seed_base = 20260619_01000),   "complete",    "complete",         :model),
    (() -> make_er(N; freq_seed_base = 20260619_02000),         "erdos_renyi", "ER",                :model),
    (() -> make_scalefree(N; gamma = 3.0, freq_seed_base = 20260619_03000), "scale_free", "scale-free", :model),
    (() -> load_empirical("animal_18_1";     label = "Bat",           dir = NETS_DIR, freq_seed_base = 20260619_04000), "animal_18_1",     "Bat",            :empirical),
    (() -> load_empirical("highschool_boys"; label = "Highschool boy",dir = NETS_DIR, freq_seed_base = 20260619_05000), "highschool_boys", "Highschool boy", :empirical),
    (() -> load_empirical("animal_6";        label = "House finch",   dir = NETS_DIR, freq_seed_base = 20260619_06000), "animal_6",        "House finch",    :empirical),
    (() -> load_empirical("surfersb";        label = "Surfer",        dir = NETS_DIR, freq_seed_base = 20260619_07000), "surfersb",        "Surfer",         :empirical),
    (() -> load_empirical("hall";            label = "Hall",          dir = NETS_DIR, freq_seed_base = 20260619_08000), "hall",            "Hall",           :empirical),
    (() -> load_empirical("jazz";            label = "Jazz",          dir = NETS_DIR, freq_seed_base = 20260619_09000), "jazz",            "Jazz",           :empirical),
]

"""
Number of sentinel nodes n = mult * floor(ln N). `mult=1` (default) is the paper's main definition; `mult=2`
is the auxiliary "doubled sentinel set size" robustness check.
"""
n_sentinels(net::Union{NetCase,NetTopology}; mult::Int = 1) = mult * floor(Int, log(net.nnodes))

"Round K to the simulation grid spacing (0.01)."
snap_to_grid(K) = round(Int, K / 0.01) * 0.01

# =============================================================================
# Sentinel-set selection
# =============================================================================
#
# Six sentinel sets are used throughout the paper:
#   S1  all N nodes (full observation)
#   S2  top-n by an EWS-specific criterion (c_i for the order-parameter
#       signals, Var_t delta_theta_i for the microscopic EWS) evaluated at
#       K = 0.9 K_c
#   S3  same criterion as S2, evaluated at K = 0.5 K_c
#   S4  top-n by degree (hub nodes)
#   S5  bottom-n by degree (peripheral nodes)
#   S6  n nodes drawn uniformly at random, redrawn independently for every
#       simulation
# S2/S3 need per-node dynamical quantities computed at specific K values
# (`node_stats_one_K`, below); S4/S5 are static topological properties
# computed once per network; S6 is redrawn per simulation. All six are
# assembled together in simulate/01_ews_sentinel_sweep.jl.

"""
    select_top(vals, n; rng)

The n indices with the largest values in `vals` (used for S2, S3, S4, S5).
Ties at the n-th-largest value are broken uniformly at random (reproducibly,
from `rng`). Exact ties are common for the degree-based sets S4/S5; for
S2/S3 (ranked by the continuous `c_i`/`Var_t delta_theta_i`) they are
essentially never tied in practice, but every caller passes its own `rng`
regardless, for consistency."""
function select_top(vals, n; rng)
    order = sortperm(vals; rev = true)
    thresh = vals[order[n]]
    above = findall(v -> v > thresh, vals)
    tied = findall(v -> v == thresh, vals)
    need = n - length(above)
    chosen = tied[randperm(rng, length(tied))[1:need]]
    sort(vcat(above, chosen))
end

"`nsets` independent sets of n indices drawn uniformly at random without replacement (S6)."
random_sets(N, n, nsets; rng) = [sort(randperm(rng, N)[1:n]) for _ in 1:nsets]

# =============================================================================
# Euler-Maruyama simulation
# =============================================================================

"""
    drift!(du, u, net, K, sbuf, cbuf, asbuf, acbuf)

In-place noiseless drift of Eq. 1. The all-to-all case uses the closed-form
mean-field sum (O(N)); all other networks use two sparse matrix-vector
products against sin(theta) and cos(theta) (O(E)). `sbuf/cbuf/asbuf/acbuf`
are pre-allocated scratch vectors reused across time steps to avoid
per-step allocation.
"""
function drift!(du, u, net::NetCase, K, sbuf, cbuf, asbuf, acbuf)
    N = length(u)
    if net.name == :complete
        Zr = 0.0; Zi = 0.0
        @inbounds @simd for j in 1:N
            Zr += cos(u[j]); Zi += sin(u[j])
        end
        invk = 1.0 / net.mean_degree
        @inbounds @simd for i in 1:N
            du[i] = net.omega[i] + K * invk * (Zi * cos(u[i]) - Zr * sin(u[i]))
        end
    else
        @inbounds @simd for i in 1:N
            sbuf[i] = sin(u[i]); cbuf[i] = cos(u[i])
        end
        mul!(asbuf, net.adjacency, sbuf)
        mul!(acbuf, net.adjacency, cbuf)
        invk = 1.0 / net.mean_degree
        @inbounds @simd for i in 1:N
            du[i] = net.omega[i] + K * invk * (cbuf[i] * asbuf[i] - sbuf[i] * acbuf[i])
        end
    end
    nothing
end

"""
    var_unbiased(sumv, sumv2, M)

Unbiased (Bessel-corrected) sample variance of M observations, given only
their running sum `sumv` and sum-of-squares `sumv2`: divides the sum of
squared deviations from the sample mean by (M - 1) rather than M. Used for
the order-parameter variance estimators Var_t R(t) and Var_t R_S(t) (both
global and per-sentinel-set), which estimate one parameter (the temporal
mean) from the M samples. Returns 0.0 if M < 2 (variance undefined/zero for
a single point).
"""
@inline function var_unbiased(sumv, sumv2, M)
    M < 2 && return 0.0
    vbar = sumv / M
    max((sumv2 - M * vbar^2) / (M - 1), 0.0)
end

"""
    resid_var(sy, syy, sxy, M)

Closed-form UNBIASED residual variance of an ordinary-least-squares fit of a
linear trend to M equally-spaced samples y_1..y_M, given only the running
sums sy=sum(y), syy=sum(y^2), sxy=sum(t*y) (t=1..M). This is the microscopic
EWS building block and lets us
compute it in a single streaming pass without storing the whole phase
trajectory.

The bracketed term below is the *population* (divide-by-M) residual
variance: Var(y) - Cov(t,y)^2/Var(t), i.e. the sum of squared OLS residuals
divided by M. Because the OLS fit estimates two parameters (slope and
intercept) from the same M points, the unbiased estimate divides that sum of
squared residuals by (M - 2) instead -- so the whole population estimate is
rescaled by M/(M - 2). Returns 0.0 if M < 3 (residual variance undefined for
fewer than 3 points once 2 degrees of freedom are used up).
"""
@inline function resid_var(sy, syy, sxy, M)
    M < 3 && return 0.0
    xbar = (M + 1) / 2
    varx = (M * M - 1) / 12.0
    ybar = sy / M
    pop_resid_var = max(syy / M - ybar^2 - (sxy / M - xbar * ybar)^2 / varx, 0.0)
    pop_resid_var * M / (M - 2)
end

"Advance the state through one K value without recording anything (used to
chain the initial condition between consecutive K values of an up-sweep)."
function advance_one_K(K, u0, net::NetCase, p::SimParams, rng)
    N = length(u0)
    dt = p.dt
    noise_amp = p.noise_sigma * sqrt(dt)
    nsteps = round(Int, p.T / dt)
    u = copy(u0)
    du = similar(u); sbuf = similar(u); cbuf = similar(u); asbuf = similar(u); acbuf = similar(u)
    @inbounds for _ in 1:nsteps
        drift!(du, u, net, K, sbuf, cbuf, asbuf, acbuf)
        for i in 1:N
            u[i] += du[i] * dt + noise_amp * randn(rng)
        end
    end
    u
end

"""
    global_one_K(K, u0, net, p, rng) -> (u_final, Rmean)

Simulate one K value and return the post-burn-in mean of the global order
parameter R(t). Used to determine the empirical K_c
(simulate/01_ews_sentinel_sweep.jl), which only ever needs <R>_t(K).
"""
function global_one_K(K, u0, net::NetCase, p::SimParams, rng)
    N = length(u0)
    dt = p.dt
    noise_amp = p.noise_sigma * sqrt(dt)
    nsteps = round(Int, p.T / dt)
    nburn = round(Int, p.burnin / dt)
    stride = round(Int, p.sample_dt / dt)
    u = copy(u0)
    du = similar(u); sbuf = similar(u); cbuf = similar(u); asbuf = similar(u); acbuf = similar(u)
    sumR = 0.0; M = 0
    @inbounds for step in 1:nsteps
        drift!(du, u, net, K, sbuf, cbuf, asbuf, acbuf)
        for i in 1:N
            u[i] += du[i] * dt + noise_amp * randn(rng)
        end
        if step > nburn && step % stride == 0
            M += 1
            zr = 0.0; zi = 0.0
            @simd for i in 1:N
                zr += cos(u[i]); zi += sin(u[i])
            end
            R = sqrt(zr * zr + zi * zi) / N
            sumR += R
        end
    end
    copy(u), sumR / M
end

"""
    sweep_global(net, Kgrid, p; seed) -> Rmean::Vector{Float64}

Chained up-sweep over `Kgrid` (one run): at each K, `global_one_K` is called
with the final state of the previous K as the initial condition (K0's
initial condition is uniform random). Returns Rmean(K).
"""
function sweep_global(net::NetCase, Kgrid, p::SimParams; seed::Int)
    rng = MersenneTwister(seed)
    nK = length(Kgrid)
    Rmean = Vector{Float64}(undef, nK)
    u = 2pi .* rand(rng, net.nnodes)
    for ik in 1:nK
        u, rm = global_one_K(Kgrid[ik], u, net, p, rng)
        Rmean[ik] = rm
    end
    Rmean
end

"""
    node_stats_one_K(K, u0, net, p, rng) -> (u_final, Rmean, vd, ci)

Simulate one K value and return, for every node i: `vd[i]` = Var_t[delta_theta_i]
(the microscopic EWS, unbiased -- see `resid_var`) and `ci[i]` =
|E_t[exp(i(Psi(t)-theta_i(t)))]| (the mean-field locking measure).
Used both for sentinel-set selection (simulate/02) and for the per-node
heatmaps.
"""
function node_stats_one_K(K, u0, net::NetCase, p::SimParams, rng)
    N = length(u0)
    dt = p.dt
    noise_amp = p.noise_sigma * sqrt(dt)
    nsteps = round(Int, p.T / dt)
    nburn = round(Int, p.burnin / dt)
    stride = round(Int, p.sample_dt / dt)
    u = copy(u0)
    du = similar(u); sbuf = similar(u); cbuf = similar(u); asbuf = similar(u); acbuf = similar(u)
    y0 = zeros(N); have_y0 = falses(N)
    sy = zeros(N); syy = zeros(N); sxy = zeros(N)
    cr = zeros(N); cs = zeros(N)
    M = 0; sumR = 0.0
    @inbounds for step in 1:nsteps
        drift!(du, u, net, K, sbuf, cbuf, asbuf, acbuf)
        for i in 1:N
            u[i] += du[i] * dt + noise_amp * randn(rng)
        end
        if step > nburn && step % stride == 0
            M += 1
            zr = 0.0; zi = 0.0
            @simd for i in 1:N
                zr += cos(u[i]); zi += sin(u[i])
            end
            sumR += sqrt(zr * zr + zi * zi) / N
            Psi = atan(zi, zr)
            @inbounds for i in 1:N
                theta = u[i]
                if !have_y0[i]
                    y0[i] = theta; have_y0[i] = true
                end
                y = theta - y0[i]
                sy[i] += y; syy[i] += y * y; sxy[i] += M * y
                a = Psi - theta
                cr[i] += cos(a); cs[i] += sin(a)
            end
        end
    end
    vd = Vector{Float64}(undef, N)
    ci = Vector{Float64}(undef, N)
    @inbounds for i in 1:N
        vd[i] = resid_var(sy[i], syy[i], sxy[i], M)
        ci[i] = sqrt(cr[i]^2 + cs[i]^2) / M
    end
    copy(u), sumR / M, vd, ci
end

"""
    record_one_K(K, u0, net, sets, p, rng)
        -> (u_final, Rmean, Rvar, global_micro, RSmean, RSvar, set_micro)

Simulate one K value and return, for a collection of node sets (sentinel
sets), each set's local order parameter R_S (`RSmean`), its UNBIASED
variance (`RSvar`), and its mean microscopic EWS (`set_micro`). Also returns the global order
parameter mean/variance and the whole-network mean microscopic EWS.
"""
function record_one_K(K, u0, net::NetCase, sets, p::SimParams, rng)
    N = length(u0)
    dt = p.dt
    noise_amp = p.noise_sigma * sqrt(dt)
    nsteps = round(Int, p.T / dt)
    nburn = round(Int, p.burnin / dt)
    stride = round(Int, p.sample_dt / dt)
    u = copy(u0)
    du = similar(u); sbuf = similar(u); cbuf = similar(u); asbuf = similar(u); acbuf = similar(u)
    nsets = length(sets)
    y0 = zeros(N); have_y0 = falses(N)
    sy = zeros(N); syy = zeros(N); sxy = zeros(N)
    sumR = 0.0; sumR2 = 0.0
    sumRS = zeros(nsets); sumRS2 = zeros(nsets)
    M = 0
    @inbounds for step in 1:nsteps
        drift!(du, u, net, K, sbuf, cbuf, asbuf, acbuf)
        for i in 1:N
            u[i] += du[i] * dt + noise_amp * randn(rng)
        end
        if step > nburn && step % stride == 0
            M += 1
            zr = 0.0; zi = 0.0
            @simd for i in 1:N
                zr += cos(u[i]); zi += sin(u[i])
            end
            R = sqrt(zr * zr + zi * zi) / N
            sumR += R; sumR2 += R * R
            @inbounds for i in 1:N
                theta = u[i]
                if !have_y0[i]
                    y0[i] = theta; have_y0[i] = true
                end
                y = theta - y0[i]
                sy[i] += y; syy[i] += y * y; sxy[i] += M * y
            end
            for (k, nodes) in enumerate(sets)
                sr = 0.0; sj = 0.0
                for j in nodes
                    sr += cos(u[j]); sj += sin(u[j])
                end
                RS = sqrt(sr * sr + sj * sj) / length(nodes)
                sumRS[k] += RS; sumRS2[k] += RS * RS
            end
        end
    end
    Rmean = sumR / M
    Rvar = var_unbiased(sumR, sumR2, M)
    vd = Vector{Float64}(undef, N)
    @inbounds for i in 1:N
        vd[i] = resid_var(sy[i], syy[i], sxy[i], M)
    end
    global_micro = mean(vd)
    RSmean = Vector{Float64}(undef, nsets)
    RSvar = similar(RSmean)
    set_micro = similar(RSmean)
    for (k, nodes) in enumerate(sets)
        RSmean[k] = sumRS[k] / M
        RSvar[k] = var_unbiased(sumRS[k], sumRS2[k], M)
        set_micro[k] = mean(@view vd[nodes])
    end
    copy(u), Rmean, Rvar, global_micro, RSmean, RSvar, set_micro
end

# =============================================================================
# Rank statistics
# =============================================================================

"Kendall's rank correlation tau between equal-length vectors x and y (tie-corrected)."
function kendall_tau(x::AbstractVector, y::AbstractVector)
    n = length(x)
    @assert n == length(y)
    nconcordant = 0; ndiscordant = 0; tiesx = 0; tiesy = 0
    @inbounds for i in 1:n-1, j in i+1:n
        dx = x[i] - x[j]; dy = y[i] - y[j]
        s = dx * dy
        if s > 0
            nconcordant += 1
        elseif s < 0
            ndiscordant += 1
        else
            dx == 0 && (tiesx += 1)
            dy == 0 && (tiesy += 1)
        end
    end
    n0 = n * (n - 1) / 2
    denom = sqrt((n0 - tiesx) * (n0 - tiesy))
    denom == 0 ? 0.0 : (nconcordant - ndiscordant) / denom
end
