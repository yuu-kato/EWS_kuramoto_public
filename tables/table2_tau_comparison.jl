# tables/table2_tau_comparison.jl
#
# Write Table (main.tex, label tab:tau_comparison): Kendall's tau, mean +-
# s.d. POOLED DIRECTLY OVER ALL 900 SIMULATIONS (9 networks x 10 frequency
# realizations x 10 noise runs), per sentinel set (S1-S6) and EWS signal
# (local order parameter R_S, collective Var_t R_S, microscopic), comparing
# the main-text sentinel-set size n=floor(ln N) against the auxiliary
# n=2*floor(ln N).
# All six sets come from the single cache/bundle_ews_<tag>[_mult2].jls
# bundle (simulate/01_ews_sentinel_sweep.jl) -- no second "hubfreq" bundle to
# merge in, unlike the earlier two-script version of this pipeline; S4 (hub)
# and S5 (periphery) are bundle columns 6 and 7, S6 (random) is column 8 --
# see simulate/01_ews_sentinel_sweep.jl's header note.
#
# For each (mult, sentinel set, signal), all 900 individual-simulation tau
# values across all 9 networks are pooled into one vector, and the "mean +-
# s.d." reported here is the mean and standard deviation of that pooled set
# -- i.e. every simulation is weighted equally regardless of which network it
# came from (the per-network view is what Figs. 1/2 already show).
#
# The optional argument selects the K_c-definition family written by
# simulate/01_ews_sentinel_sweep.jl: Kc20 (default; the paper's threshold
# fraction 0.2), or Kc10 / Kc30 (fractions 0.1 / 0.3, robustness checks of
# the K_c definition). Kendall's tau is always evaluated over K_0 <= K <= K_c
# of that family, and S2/S3 are selected at 0.9/0.5 of that K_c.
#
# Run (after simulate/01_ews_sentinel_sweep.jl with the matching fraction,
# which writes both bundle_ews_<tag>.jls and bundle_ews_<tag>_mult2.jls in
# one run):
#   julia --project=. tables/table2_tau_comparison.jl          # Kc20 (paper)
#   julia --project=. tables/table2_tau_comparison.jl Kc10
#   julia --project=. tables/table2_tau_comparison.jl Kc30
#
# Writes ../draft/table2_tau_comparison[_Kc10|_Kc30].tex (same directory as
# main.tex) -- just the table itself (\label + ruledtabular + tabular), no
# \begin{table*}/\caption wrapper; \input this file inside that float.

include(joinpath(@__DIR__, "..", "core", "kuramoto_core.jl"))
using Serialization, Statistics, Printf

const HERE  = @__DIR__
const CACHE = joinpath(HERE, "..", "cache")
const OUTDIR = normpath(joinpath(HERE, "..", "draft")); mkpath(OUTDIR)
const SET_LABELS = ["S1", "S2", "S3", "S4", "S5", "S6"]
# (bundle field, signal symbol, column header)
const SIGNALS = [(:tau_RS, :RS, "\$\\langle R_S\\rangle_t\$"), (:tau_RSvar, :RSvar, "\$\\mathrm{Var}(R_S)\$"),
                  (:tau_micro, :micro, "\$V_S\$")]
const MULTS = [1, 2]
const KCTAG = length(ARGS) >= 1 ? ARGS[1] : "Kc20"   # bundle family to read: Kc10, Kc20, Kc30
const KCSUF = KCTAG == "Kc20" ? "" : "_" * KCTAG

"""Bundle column (1..8, S1/S2-1/S2-2/S3-1/S3-2/S4/S5/S6 -- see
simulate/01_ews_sentinel_sweep.jl's header note) for sentinel set `newset`
(1..6, paper numbering) and `signal`
-- identical mapping to simulate/02_tipmoc_export_csv.jl's `newset_column`
and figures/fig_ews_sentinel.jl's `set_column`."""
function newset_column(newset, signal)
    newset == 1 && return 1
    newset == 2 && return signal == :micro ? 3 : 2
    newset == 3 && return signal == :micro ? 5 : 4
    newset == 4 && return 6   # S4: hub
    newset == 5 && return 7   # S5: periphery
    newset == 6 && return 8   # S6: random
    error("newset must be 1..6")
end

"All 900 individual-simulation tau values (9 networks x 10 frequency
realizations x 10 noise runs), pooled across networks, for sentinel set
`newset` (1..6) and `signal`, from bundle `B`."
function pooled_tau_values(B, fld, newset, signal)
    col = newset_column(newset, signal)
    reduce(vcat, [getfield(rr, fld)[col] for r in 1:B.nfreq for rr in B.results_allfreq[i][r]]
                 for i in eachindex(B.meta))
end

"Mean and s.d. of tau POOLED DIRECTLY OVER ALL 900 SIMULATIONS, for sentinel
set `newset` (1..6) and `signal`, from bundle `B`."
function pooled_mean_sd(B, fld, newset, signal)
    vals = pooled_tau_values(B, fld, newset, signal)
    mean(vals), std(vals)
end

bundle(mult) = deserialize(joinpath(CACHE, mult == 1 ? "bundle_ews_$(KCTAG).jls" : "bundle_ews_$(KCTAG)_mult$(mult).jls"))
B = Dict(mult => bundle(mult) for mult in MULTS)

"Format a mean value, wrapping negative numbers in math mode so LaTeX
renders a proper minus sign instead of a plain hyphen."
fmt_mean(mu) = mu < 0 ? "\$" * @sprintf("%.2f", mu) * "\$" : @sprintf("%.2f", mu)

"Cell for sentinel set `newset` at (mult, fld, signal)."
function cell(mult, newset, fld, signal)
    mu, sd = pooled_mean_sd(B[mult], fld, newset, signal)
    @sprintf("%s\\,(%.2f)", fmt_mean(mu), sd)
end

io = IOBuffer()
println(io, "\\label{tab:tau_comparison}")
println(io, "\\begin{ruledtabular}")
println(io, "\\begin{tabular}{l", repeat("c", length(SIGNALS) * length(MULTS)), "}")
print(io, " & ")
println(io, join(["\\multicolumn{$(length(SIGNALS))}{c}{\$n=$(mult == 1 ? "\\lfloor\\ln N\\rfloor" : "2\\lfloor\\ln N\\rfloor")\$}" for mult in MULTS], " & "), " \\\\")
print(io, " & ")
println(io, join(fill(join([hdr for (_, _, hdr) in SIGNALS], " & "), length(MULTS)), " & "), " \\\\")
println(io, "\\colrule")
for (newset, label) in enumerate(SET_LABELS)
    cells = String[cell(mult, newset, fld, signal) for mult in MULTS for (fld, signal, _) in SIGNALS]
    println(io, label, " & ", join(cells, " & "), " \\\\")
end
println(io, "\\end{tabular}")
println(io, "\\end{ruledtabular}")
write(joinpath(OUTDIR, "table2_tau_comparison$(KCSUF).tex"), String(take!(io)))
println("wrote table2_tau_comparison$(KCSUF).tex to ", OUTDIR)
