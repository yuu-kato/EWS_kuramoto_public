# figures/fig_ews_sentinel_appendix.jl
#
# Appendix counterpart to figures/fig_ews_sentinel.jl (main.tex Fig. 2): the
# same three-EWS-vs-K/K_c panels, but for the six networks NOT shown in
# the main text -- ER, Highschool boy, House finch, Surfer, Hall, Jazz,
# anchored to K_c (cache/bundle_ews_Kc20.jls). Sentinel sets S2/S3 use a
# signal-dependent node selection, S4/S5 are the static hub/periphery sets,
# and S6 is n nodes chosen uniformly at random -- a single independent draw
# per simulation, the same node set used for all three signals -- see
# figures/fig_ews_sentinel.jl's header note for the full explanation and
# `newset_column` for the exact per-signal column mapping (no new simulation
# needed). Illustrative figure: shows a single representative run (run 2) of
# the first frequency realization's NRUN_PER_FREQ=10 noise runs, raw (not
# averaged/smoothed). S1 (all nodes) is deliberately omitted here -- it's the
# whole-network reference already shown in fig_ews_all_kc_appendix.jl, so
# repeating it in every sentinel-set panel would be redundant.
# Legend (S2..S6 only) sits inside the first panel.
# Panel labels (A1)-(A3) .. (F1)-(F3) (one letter per network row, one
# number per EWS column) are placed outside each panel's top-left corner.
# Only the K_c reference line is drawn (dotted); the 0.5*K_c reference is
# omitted, as in fig_ews_sentinel.jl.
#
# Run (after simulate/01_ews_sentinel_sweep.jl):
#   julia --project=. figures/fig_ews_sentinel_appendix.jl
#
# Writes ../draft/ews_sentinel_appendix.pdf.

include(joinpath(@__DIR__, "..", "core", "kuramoto_core.jl"))   # defines SimParams, needed to deserialize bundle_ews_Kc20.jls
include(joinpath(@__DIR__, "..", "core", "figure_style.jl")); using .FigStyle
using CairoMakie, Serialization, Statistics, Printf
CairoMakie.update_theme!(Axis = (; xgridvisible = false, ygridvisible = false))

const HERE  = @__DIR__
const CACHE = joinpath(HERE, "..", "cache")
const OUTDIR = normpath(joinpath(HERE, "..", "draft")); mkpath(OUTDIR)
const B = deserialize(joinpath(CACHE, "bundle_ews_Kc20.jls"))
const KEYS = ["erdos_renyi", "highschool_boys", "animal_6", "surfersb", "hall", "jazz"]
const REALIZATION = 1
const RUN = 2   # single representative run, shown raw (see header note)
const ROWS = [findfirst(m -> String(m.key) == k, B.meta) for k in KEYS]

# ---- style knobs ----
const NNET = length(ROWS)
const ROWHEIGHT = 420
const FIGSIZE = (2200, 80 + ROWHEIGHT * NNET)
const PANELLABELSIZE = 40
const TICKLABELSIZE = 32
const SET_LABELS = ["S1", "S2", "S3", "S4", "S5", "S6"]   # legend needs only the short set names
# S1's slot (black) is kept reserved but unused, so S2..S6 keep the same
# colour/marker identity as in any other figure that does include S1:
const PAL = [RGBf(0, 0, 0), RGBf(0.0, 0.45, 0.70), RGBf(0.90, 0.45, 0.0), RGBf(0.0, 0.60, 0.50),
             RGBf(0.58, 0.20, 0.60), RGBf(0.85, 0.10, 0.10)]
const MARKERS = [:circle, :rect, :utriangle, :diamond, :star5, :xcross]
const OTHER_LW = 4.0; const BAND_ALPHA = 0.13
const MARKERSIZE = 30; const N_MARKS = 9
const KC_LW = 3.0

"""Bundle column (1..8, S1/S2-1/S2-2/S3-1/S3-2/S4/S5/S6) for sentinel set
`newset` in 2:6 and `signal` -- see figures/fig_ews_sentinel.jl's identical
helper for the full explanation."""
function newset_column(newset, signal)
    newset == 2 && return signal == :micro ? 3 : 2
    newset == 3 && return signal == :micro ? 5 : 4
    newset == 4 && return 6   # S4: hub
    newset == 5 && return 7   # S5: periphery
    newset == 6 && return 8   # S6: random
    error("newset must be 2..6 (S1 is handled separately, see fig_ews_all_kc_appendix.jl)")
end

"Column c of `field` for run `rr`, for sentinel set `newset` (2..6) and
`signal`. The macroscopic field (RSvar) is rescaled by the set size n to
make it intensive (Var_t R_S ~ 1/n), per main.tex Sec. 'Early warning
signals'."
function set_column(rr, field, newset, signal)
    col = newset_column(newset, signal)
    values = getfield(rr, field)[:, col]
    field === :RSvar ? values .* length(rr.sets[col]) : values
end

"`Krel` is K/K_c; K_c itself is always at 1 once rescaled. `rr` is the
single representative run's record (raw trajectory, not averaged across
runs -- see header note)."
function ews_panel!(ax, rr, Krel, field, signal; show_legend = false)
    for newset in 2:6                     # S1 (all nodes) omitted -- see fig_ews_all_kc_appendix.jl instead
        vals = set_column(rr, field, newset, signal)
        lines!(ax, Krel, vals; color = PAL[newset], linewidth = OTHER_LW, label = SET_LABELS[newset])
    end
    vlines!(ax, [1.0]; color = :black, linestyle = :dot, linewidth = KC_LW)
    xlims!(ax, 0, 2.0)
    show_legend && axislegend(ax; position = :rt, framevisible = false, labelsize = FigStyle.LEGENDSIZE,
                               nbanks = 1, merge = true, patchsize = (70, 3))
end

"Row letter for panel labels: row 1 -> A, row 2 -> B, ..."
row_letter(row) = string(Char('A' + row - 1))

# (field, signal, ylabel)
const COLUMNS = [
    (:RS,    :RS,    L"\langle R_S \rangle_t"),
    (:RSvar, :RSvar, L"n \,\, \mathrm{Var}(R_S)"),
    (:SM,    :micro, L"V_S"),
]

fig = Figure(size = FIGSIZE, figure_padding = (20, 50, 10, 10))   # extra right padding: last x-tick label ("2.0") sits right at the panel edge

for (row, i) in enumerate(ROWS)
    m = B.meta[i]; rm = B.realization_meta[i][REALIZATION]
    rr = B.results_allfreq[i][REALIZATION][RUN]
    Krel = rr.Kv ./ rm.Kc

    Label(fig[row, 0], m.label; rotation = pi/2, tellheight = false, fontsize = ROWLABELSIZE)

    for (col, (field, signal, ylab)) in enumerate(COLUMNS)
        ax = FigStyle.styled_axis(fig[row, col];
            xlabel = row == NNET ? L"K/K_\mathrm{c}" : "", ylabel = ylab,
            xticklabelsize = TICKLABELSIZE, xticklabelsvisible = row == length(ROWS), yticklabelsize = TICKLABELSIZE)
        Label(fig[row, col, TopLeft()], "($(row_letter(row))$(col))";
              fontsize = PANELLABELSIZE, font = :bold, padding = (0, 14, 12, 0), halign = :left)
        ews_panel!(ax, rr, Krel, field, signal; show_legend = (row == 1 && col == 3))
        # shared range across all 9 networks -- see fig_ews_sentinel.jl's identical comment
        col == 1 && ylims!(ax, 0.25, 1.05)
        col == 2 && ylims!(ax, -0.02, 0.42)
        col == 3 && ylims!(ax, -0.5, 16.0)
    end
end

colsize!(fig.layout, 0, Fixed(30))
for col in 1:length(COLUMNS)
    colsize!(fig.layout, col, Relative(1 / length(COLUMNS)))
end
colgap!(fig.layout, 60)

save(joinpath(OUTDIR, "ews_sentinel_appendix.pdf"), fig)
println("wrote ews_sentinel_appendix.pdf to ", OUTDIR)
