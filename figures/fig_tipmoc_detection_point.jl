# figures/fig_tipmoc_detection_point.jl
#
# Fig. (main.tex, label fig:tipmoc_detection_point): distribution of the
# TIPMOC detection point K_det / K_c, for each sentinel set S1-S6 (rows) and
# each EWS signal (columns: local order parameter R_S, collective Var_t R_S,
# microscopic V_S), pooled over all 9 networks x 10 frequency realizations x
# 10 noise runs = 900 trajectories per panel, at the main-text sentinel-set
# size n = floor(ln N). K_c is the p=20% threshold-crossing definition (Kc20).
#
# Unlike Table 3 (tables/table3_tipmoc_detection.jl), which only asks whether
# TIPMOC fires BEFORE K_c, this figure uses the detector output for the whole
# stored sweep K_0..2K_c (simulate/02_tipmoc_export_csv.jl run with `full`,
# then simulate/02_tipmoc_detect.py with `full`), so a trajectory can be
# detected before K_c (blue), detected only at or after K_c (orange -- the
# detector reacting to the transition itself, not an early warning), or never
# detected within 2K_c. The annotation in each panel gives this split in
# whole percent ('< 1%' for nonzero shares below 1%); the blue share is the
# detection rate reported in Table 3. The last line is the median of K_det
# over the detections before K_c, in units of K_c, i.e. how early TIPMOC
# typically fires when it does give an early warning. The scan
# stops at the first detection, so K_det is the first K at which TIPMOC
# declares a transition (column u_det of the results CSV). Each panel is a
# histogram of K_det / K_c with bin width 0.05, normalised by the 900 runs so
# that every panel shares one y axis (fraction of runs). Dotted line: K_c.
#
# Layout matches Figs. 2/3 (fig_ews_all_kc.jl, fig_ews_sentinel.jl)
# transposed: 6 rows (S1..S6, top to bottom) x 3 columns (R_S, Var_t R_S,
# V_S, left to right). Panel labels (A1)-(A3), ..., (F1)-(F3). Legend sits
# above the panels.
#
# Run (after simulate/02_tipmoc_detect.py 1 full):
#   julia --project=. figures/fig_tipmoc_detection_point.jl
#
# Writes ../draft/tipmoc_detection_point.pdf.

include(joinpath(@__DIR__, "..", "core", "figure_style.jl")); using .FigStyle
using CairoMakie, DelimitedFiles, Statistics, Printf
CairoMakie.update_theme!(Axis = (; xgridvisible = false, ygridvisible = false))

const HERE  = @__DIR__
const CACHE = joinpath(HERE, "..", "cache")
const OUTDIR = normpath(joinpath(HERE, "..", "draft")); mkpath(OUTDIR)
const RESULTS = joinpath(CACHE, "tipmoc_results_Kc20_to2Kc.csv")   # mult=1, K_0..2K_c window

# ---- style knobs ----
const FIGSIZE = (2400, 3000)
const ROWLABELSIZE = 50
const PANELLABELSIZE = 50
const TICKLABELSIZE = 36
const ANNOTSIZE = 40
const SET_LABELS = ["S1", "S2", "S3", "S4", "S5", "S6"]
const BINS = 0.0:0.05:2.0      # in units of K/K_c
const YMAX = 0.28              # fraction of runs; shared by every panel
const C_EARLY = RGBf(0.10, 0.30, 0.55)   # detected before K_c
const C_LATE  = RGBf(0.85, 0.45, 0.10)   # detected at or after K_c
const KC_LW = 3.0

# ---- data: simulate/02_tipmoc_detect.py results CSV, read as in tables/table3 ----
raw, hdr = readdlm(RESULTS, ','; header = true)
hdr = vec(hdr)
col(name) = raw[:, findfirst(==(name), hdr)]
scopes   = String.(col("scope")); signals = String.(col("signal"))
detected = String.(col("detected")) .== "True"
u_det    = Float64.(replace(col("u_det"), "" => NaN))
Kc       = Float64.(col("Kc"))
ratio    = u_det ./ Kc     # K_det / K_c (NaN where not detected)

"Percentage label: whole percent, except '< 1%' for nonzero shares below 1%."
pct(x) = x == 0 ? "0%" : x < 1 ? "< 1%" : @sprintf("%.0f%%", x)

"Histogram panel for sentinel set `set` and results-CSV signal name
`signal` ('R_S', 'VarR', 'micro'): detections before K_c in blue, at or
after K_c in orange, both normalised by the total number of runs in the
cell; annotation gives the before/after/never percentages."
function detection_panel!(ax, set, signal)
    cell  = (scopes .== set) .& (signals .== signal)
    ntot  = sum(cell)
    r     = ratio[cell .& detected]
    early = r[r .< 1 - 1e-9]; late = r[r .>= 1 - 1e-9]
    isempty(early) || hist!(ax, early; bins = BINS, weights = fill(1 / ntot, length(early)), color = C_EARLY, strokewidth = 0)
    isempty(late)  || hist!(ax, late;  bins = BINS, weights = fill(1 / ntot, length(late)),  color = C_LATE,  strokewidth = 0)
    vlines!(ax, [1.0]; color = :black, linestyle = :dot, linewidth = KC_LW)
    isempty(early) || vlines!(ax, [median(early)]; color = C_EARLY, linestyle = :dash, linewidth = KC_LW+1)
    xlims!(ax, 0, 2.0); ylims!(ax, 0, YMAX)
    med_line = isempty(early) ? "" : @sprintf("(median %.2f)\n", median(early))   # omitted when nothing is detected before K_c
    text!(ax, 0.97, 0.95; space = :relative, align = (:right, :top), fontsize = ANNOTSIZE,
          text = string("before ", pct(100 * length(early) / ntot), "\n",
                        med_line,
                        "after ",  pct(100 * length(late) / ntot), "\n",
                        "never ",  pct(100 * (ntot - length(r)) / ntot)))
end

row_letter(row) = string(Char('A' + row - 1))

# (results-CSV signal name, column title)
const COLS_SPEC = [
    ("R_S",   L"\langle R_S \rangle_t"),
    ("VarR",  L"\mathrm{Var}(R_S)"),
    ("micro", L"V_S"),
]

fig = Figure(size = FIGSIZE, figure_padding = (20, 30, 10, 10))

Legend(fig[0, 1:length(COLS_SPEC)],
       [PolyElement(color = C_EARLY), PolyElement(color = C_LATE)],
       [L"\text{detected before } K_\mathrm{c}", L"\text{detected at or after } K_\mathrm{c}"];
       orientation = :horizontal, framevisible = false, labelsize = FigStyle.LEGENDSIZE + 5)

for (row, set) in enumerate(SET_LABELS)
    Label(fig[row, 0], set; rotation = pi/2, tellheight = false, fontsize = ROWLABELSIZE)

    for (col, (signal, title)) in enumerate(COLS_SPEC)
        ax = FigStyle.styled_axis(fig[row, col];
            title = row == 1 ? title : "",
            titlesize = 50, 
            xlabel = row == length(SET_LABELS) ? L"K_\mathrm{det}/K_\mathrm{c}" : "",
            xlabelsize = 50, 
            ylabel = col == 1 ? "fraction of runs" : "",
            ylabelsize = 45, 
            xticklabelsvisible = row == length(SET_LABELS), yticklabelsvisible = col == 1,
            xticklabelsize = TICKLABELSIZE, yticklabelsize = TICKLABELSIZE,
            xticks = ([0, 0.5, 1.0, 1.5, 2.0], ["0", "0.5", "1", "1.5", "2"]),
            yticks = collect(0.0:0.1:0.20))
        Label(fig[row, col, TopLeft()], "($(row_letter(row))$(col))";
              fontsize = PANELLABELSIZE, font = :bold, padding = (80, -80, 5, 0), halign = :left)
        detection_panel!(ax, set, signal)
    end
end

colsize!(fig.layout, 0, Fixed(40))
for col in 1:length(COLS_SPEC)
    colsize!(fig.layout, col, Relative(1 / length(COLS_SPEC)))
end
colgap!(fig.layout, 50)
rowgap!(fig.layout, 30)

save(joinpath(OUTDIR, "tipmoc_detection_point.pdf"), fig)
println("wrote tipmoc_detection_point.pdf to ", OUTDIR)
