# tables/table1_networks_and_Kc.jl
#
# Write Table I of main.tex (label tab:kc): the nine networks with N, E,
# <k>, n = floor(ln N), and the order-parameter-threshold
# characterization of the critical coupling, Kc20 (mean +- s.d. over
# the NFREQ=10 independent frequency realizations, each with its own value
# -- see core/kuramoto_core.jl's `Kc_threshold` /
# simulate/01_ews_sentinel_sweep.jl). 
#
# Run (after simulate/01_ews_sentinel_sweep.jl):
#   julia --project=. tables/table1_networks_and_Kc.jl
#
# Writes ../draft/table1_networks_and_Kc.tex (same directory as main.tex) --
# just the table itself (\label + ruledtabular + tabular), no
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

m, hdr = readdlm(joinpath(CACHE, "kc_summary.tsv"), '\t', header = true)
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
    Kc20_mean, Kc20_sd = Float64(m[i, findfirst(==("Kc20_mean"), hdr)]), Float64(m[i, findfirst(==("Kc20_sd"), hdr)])
    if i > 1 && kind != lastkind
        println(io, "\\colrule")
    end
    global lastkind = kind
    @printf(io, "%-14s & %4d & %6d & %6.2f & %d & %.2f\\,(%.2f) \\\\\n",
        label, N, E, kbar, floor(Int, log(N)), Kc20_mean, Kc20_sd)
end
println(io, "\\end{tabular}")
println(io, "\\end{ruledtabular}")
write(joinpath(OUTDIR, "table1_networks_and_Kc_v7.tex"), String(take!(io)))
println("wrote table1_networks_and_Kc_v7.tex to ", OUTDIR)
