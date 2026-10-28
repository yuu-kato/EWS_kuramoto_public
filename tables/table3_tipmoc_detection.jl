# tables/table3_tipmoc_detection.jl
#
# Write Table (main.tex, label tab:tipmoc_detection): TIPMOC detection rate
# (%) per sentinel set (S1-S6) and EWS signal, pooled over all 9 networks and
# all 100 (10 frequency realizations x 10 noise runs) trajectories per cell,
# anchored to the K_c definition selected by the optional argument (default Kc20), for both
# sentinel-set sizes n=floor(ln N) and n=2*floor(ln N) side by side (same
# both-n layout as tables/table2_tau_comparison.jl).
#
# TIPMOC detection itself runs in the original Python/SciPy implementation
# (simulate/02_tipmoc_detect.py, using tipmoc_lib/delta_aicc.py and
# fit_powerlaw.py directly -- not a Julia re-implementation of them; see that
# script's header note for why), reading cache/tipmoc_results_<tag>[_mult<mult>][_to2Kc].csv
# (simulate/02_tipmoc_export_csv.jl, from cache/bundle_ews_<tag>[_mult<mult>].jls).
# This table reads that Python script's per-(network,scope,signal,run) output,
# cache/tipmoc_results_<tag>[_mult<mult>][_to2Kc].csv, and pools it into the per-cell
# detection rate reported here.
#
# Run (after simulate/02_tipmoc_export_csv.jl and simulate/02_tipmoc_detect.py,
# both mult=1 and mult=2). For each sentinel size the script reads the
# K0..2Kc detector output (tipmoc_results_<tag>[_mult2]_to2Kc.csv) if it
# exists, else the K < K_c output; only detections with u_det < K_c count, so
# both give the same statistic. The optional argument selects the K_c
# definition (Kc10/Kc20/Kc30; default Kc20):
#   julia --project=. tables/table3_tipmoc_detection.jl
#   julia --project=. tables/table3_tipmoc_detection.jl Kc10
#
# Writes ../draft/table3_tipmoc_detection[_Kc10|_Kc30].tex (same directory as main.tex)
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
const KCTAG = length(ARGS) >= 1 ? ARGS[1] : "Kc20"                      # K_c definition: Kc10, Kc20, Kc30
const KCSUF = KCTAG == "Kc20" ? "" : "_" * KCTAG

"Cross-network mean detection rate (%) for `scope`/`signal`, from a
simulate/02_tipmoc_detect.py results CSV. A run counts as detected only if
u_det < K_c, i.e. before the transition; for the plain K < K_c input this is
identical to `detected`, for the K0..2Kc input it drops detections at or
after K_c. For each network, pool over its 100 runs into a rate, then average
those 9 per-network rates -- equals the pooled-900 rate since every network
contributes exactly 100 runs."
function detection_rates(mult)
    tag = mult == 1 ? "" : "_mult$(mult)"
    f_full = joinpath(CACHE, "tipmoc_results_$(KCTAG)$(tag)_to2Kc.csv")   # K0..2Kc run, if it was done
    f = isfile(f_full) ? f_full : joinpath(CACHE, "tipmoc_results_$(KCTAG)$(tag).csv")
    raw, hdr = readdlm(f, ','; header = true)
    hdr = vec(hdr)
    col(name) = raw[:, findfirst(==(name), hdr)]
    nets = String.(col("net")); scopes = String.(col("scope")); signals = String.(col("signal"))
    detected = String.(col("detected")) .== "True"

    u_det = Float64.(replace(col("u_det"), "" => NaN))
    Kc    = Float64.(col("Kc"))
    before_Kc = u_det .< Kc .- 1e-9
    early = detected .& before_Kc
    # late  = detected .& .!before_Kc

    rates = Dict{Tuple{String,String},Float64}()
    for (signal, _) in SIGNALS, label in SET_LABELS
        cell = (scopes .== label) .& (signals .== signal)
        net_rates = [100 * mean(early[cell .& (nets .== net)]) for net in unique(nets)]
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
write(joinpath(OUTDIR, "table3_tipmoc_detection$(KCSUF).tex"), String(take!(io)))
println("wrote table3_tipmoc_detection$(KCSUF).tex to ", OUTDIR)
