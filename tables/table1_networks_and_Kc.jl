# tables/table1_networks_and_Kc.jl
#
# Write Table I of main.tex (label tab:kc): the nine networks with N, E,
# <k>, n = floor(ln N), and the order-parameter-threshold
# characterization of the critical coupling, K_c (mean +- s.d. over
# the NFREQ=10 independent frequency realizations, each with its own value
# -- see core/kuramoto_core.jl's `Kc_threshold` /
# simulate/01_ews_sentinel_sweep.jl).
#
# The optional argument selects the K_c-definition family written by
# simulate/01_ews_sentinel_sweep.jl: Kc20 (default; the paper's threshold
# fraction 0.2), or Kc10 / Kc30 (fractions 0.1 / 0.3, robustness checks of
# the K_c definition). It picks cache/kc_summary[_Kc10|_Kc30].tsv as input
# and tags the output file name accordingly.
#
# Run (after simulate/01_ews_sentinel_sweep.jl with the matching fraction):
#   julia --project=. tables/table1_networks_and_Kc.jl          # Kc20 (paper)
#   julia --project=. tables/table1_networks_and_Kc.jl Kc10
#   julia --project=. tables/table1_networks_and_Kc.jl Kc30
#
# Writes ../draft/table1_networks_and_Kc[_Kc10|_Kc30].tex (same directory as
# main.tex) -- just the table itself (\label + ruledtabular + tabular), no
# \begin{table*}/\caption wrapper. The caption is
# written by hand in main.tex; \input this file inside that float, e.g.:
#   \begin{table*}[t]
#   \centering
#   \caption{...}
#   \input{table1_networks_and_Kc.tex}
#   \end{table*}

using DelimitedFiles, Printf

const HERE  = @__DIR__
const CACHE = joinpath(HERE, "..", "cache")
const OUTDIR = normpath(joinpath(HERE, "..", "draft")); mkpath(OUTDIR)
const KCTAG = length(ARGS) >= 1 ? ARGS[1] : "Kc20"   # Kc10, Kc20, Kc30
const KCSUF = KCTAG == "Kc20" ? "" : "_" * KCTAG

m, hdr = readdlm(joinpath(CACHE, "kc_summary$(KCSUF).tsv"), '\t', header = true)
hdr = vec(hdr)
col(name) = m[:, findfirst(==(name), hdr)]

io = IOBuffer()
println(io, "\\label{tab:kc}")
println(io, "\\begin{ruledtabular}")
println(io, "\\begin{tabular}{lrrrrr}")
println(io, "Network & \$N\$ & \$E\$ & \$\\langle k\\rangle\$ & \$\\lfloor\\ln N\\rfloor\$ & \$K_\\mathrm{c}\$ \\\\")
println(io, "\\colrule")
lastkind = ""
for i in 1:size(m, 1)
    label, kind = string(m[i, findfirst(==("label"), hdr)]), string(m[i, findfirst(==("kind"), hdr)])
    N, E = Int(m[i, findfirst(==("N"), hdr)]), Int(m[i, findfirst(==("E"), hdr)])
    kbar = Float64(m[i, findfirst(==("kbar"), hdr)])
    Kc_mean, Kc_sd = Float64(m[i, findfirst(==("Kc_mean"), hdr)]), Float64(m[i, findfirst(==("Kc_sd"), hdr)])
    if i > 1 && kind != lastkind
        println(io, "\\colrule")
    end
    global lastkind = kind
    @printf(io, "%-14s & %4d & %6d & %6.2f & %d & %.2f\\,(%.2f) \\\\\n",
        label, N, E, kbar, floor(Int, log(N)), Kc_mean, Kc_sd)
end
println(io, "\\end{tabular}")
println(io, "\\end{ruledtabular}")
write(joinpath(OUTDIR, "table1_networks_and_Kc$(KCSUF).tex"), String(take!(io)))
println("wrote table1_networks_and_Kc$(KCSUF).tex to ", OUTDIR)
