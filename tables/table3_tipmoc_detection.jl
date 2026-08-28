# tables/table3_tipmoc_detection.jl
#
# Write Table (main.tex, label tab:tipmoc_detection): TIPMOC detection rate
# (%) per sentinel set (S1-S6) and EWS signal, pooled over all 9 networks and
# all 100 (10 frequency realizations x 10 noise runs) trajectories per cell,
# anchored to K_c^{20%} (the sole anchor used throughout the paper), for both
# sentinel-set sizes n=floor(ln N) and n=2*floor(ln N) side by side (same
# both-n layout as tables/table2_tau_comparison.jl).
#
# TIPMOC detection itself runs in the original Python/SciPy implementation
# (simulate/02_tipmoc_detect.py, using tipmoc_lib/delta_aicc.py and
# fit_powerlaw.py directly -- not a Julia re-implementation of them; see that
# script's header note for why), reading cache/tipmoc_input_Kc20[_mult<mult>].csv
# (simulate/02_tipmoc_export_csv.jl, from cache/bundle_ews_Kc20[_mult<mult>].jls).
# This table reads that Python script's per-(network,scope,signal,run) output,
# cache/tipmoc_results_Kc20[_mult<mult>].csv, and pools it into the per-cell
# detection rate reported here.
#
# Run (after simulate/02_tipmoc_export_csv.jl and simulate/02_tipmoc_detect.py,
# both mult=1 and mult=2):
#   julia --project=. tables/table3_tipmoc_detection.jl
#
# Writes ../draft/table3_tipmoc_detection.tex (same directory as main.tex)
# -- just the table itself (\label + ruledtabular + tabular), no
# \begin{table}/\caption wrapper. The caption is
# written by hand in main.tex; \input this file inside that float, e.g.:
#   \begin{table}[t]
#   \centering
#   \caption{...}
#   \input{table3_tipmoc_detection.tex}
#   \end{table}

using DelimitedFiles, Statistics, Printf

const HERE  = @__DIR__
const CACHE = joinpath(HERE, "..", "cache")
const OUTDIR = normpath(joinpath(HERE, "..", "draft")); mkpath(OUTDIR)
const SET_LABELS = ["S1", "S2", "S3", "S4", "S5", "S6"]
const SIGNALS = [("R_S", "\$\\langle R_S\\rangle_t\$"), ("VarR", "\$\\mathrm{Var}(R_S)\$"), ("micro", "\$V_S\$")]
const MULTS = [1, 2]

"Cross-network mean detection rate (%) for `scope`/`signal`, from a
simulate/02_tipmoc_detect.py results CSV: for each network, pool `detected`
over that network's 100 runs into a rate, then average those 9 per-network
rates -- equals the pooled-900 rate since every network contributes exactly
100 runs."
function detection_rates(mult)
    tag = mult == 1 ? "" : "_mult$(mult)"
    raw, hdr = readdlm(joinpath(CACHE, "tipmoc_results_Kc20$(tag).csv"), ','; header = true)
    hdr = vec(hdr)
    col(name) = raw[:, findfirst(==(name), hdr)]
    nets = String.(col("net")); scopes = String.(col("scope")); signals = String.(col("signal"))
    detected = String.(col("detected")) .== "True"

    rates = Dict{Tuple{String,String},Float64}()
    for (signal, _) in SIGNALS, label in SET_LABELS
        cell = (scopes .== label) .& (signals .== signal)
        net_rates = [100 * mean(detected[cell .& (nets .== net)]) for net in unique(nets)]
        rates[(label, signal)] = mean(net_rates)
    end
    rates
end

R = Dict(mult => detection_rates(mult) for mult in MULTS)

io = IOBuffer()
println(io, "\\label{tab:tipmoc_detection}")
println(io, "\\begin{ruledtabular}")
println(io, "\\begin{tabular}{l", repeat("c", length(SIGNALS) * length(MULTS)), "}")
print(io, " & ")
println(io, join(["\\multicolumn{$(length(SIGNALS))}{c}{\$n=$(mult == 1 ? "\\lfloor\\ln N\\rfloor" : "2\\lfloor\\ln N\\rfloor")\$}" for mult in MULTS], " & "), " \\\\")
print(io, " & ")
println(io, join(fill(join([hdr for (_, hdr) in SIGNALS], " & "), length(MULTS)), " & "), " \\\\")
println(io, "\\colrule")
for label in SET_LABELS
    cells = String[@sprintf("%.0f", R[mult][(label, signal)]) for mult in MULTS for (signal, _) in SIGNALS]
    println(io, label, " & ", join(cells, " & "), " \\\\")
end
println(io, "\\end{tabular}")
println(io, "\\end{ruledtabular}")
write(joinpath(OUTDIR, "table3_tipmoc_detection.tex"), String(take!(io)))
println("wrote table3_tipmoc_detection.tex to ", OUTDIR)
