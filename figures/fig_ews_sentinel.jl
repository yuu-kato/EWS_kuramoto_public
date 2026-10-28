# figures/fig_ews_sentinel.jl
#
# Fig. 2 (main.tex, label fig:ews): the three EWS (local order parameter R_S,
# macroscopic n*Var_t R_S, microscopic mean_{i in S} Var_t delta_theta_i)
# versus K/K_c, for the three representative networks shown in the main
# text -- Complete, Scale-free, Bat (the other six networks are shown in
# Appendix figures/fig_ews_sentinel_appendix.jl) -- and sentinel sets S2--S6
# (colour/marker within each panel), anchored to K_c (the p=20%
# threshold-crossing definition, cache/bundle_ews_Kc20.jls; see main.tex Eq.
# def_kcp). S2 and S3 use a signal-dependent node selection (main.tex Sec.
# "Selection of sentinel nodes"): the top-n nodes by locking measure c_i for
# R_S and Var_t R_S, or by microscopic variance Var_t delta_theta_i for the
# microscopic EWS, evaluated at K=0.9K_c (S2) or K=0.5K_c (S3). S4 (hub,
# top-n degree) and S5 (periphery, bottom-n degree) are static topological
# sets, not signal-dependent. S6 is n nodes chosen uniformly at random -- a
# single independent draw per simulation, the same node set used for all
# three signals. All six live in one bundle (simulate/01_ews_sentinel_sweep.jl,
# which tracks S1-S6 together from the same trajectories -- unlike an
# earlier two-script/two-bundle version of this pipeline), so no new
# simulation is needed here -- see `newset_column` below for the exact
# per-signal column mapping. Illustrative figure: shows a single
# representative run (run 2) of the first frequency realization's
# NRUN_PER_FREQ=10 noise runs, raw (not averaged/smoothed), so the reader
# sees what one actual trajectory looks like rather than a run-averaged
# curve. S1 (all nodes) is deliberately omitted here -- it's the
# whole-network reference already shown in fig_ews_all_kc.jl, so repeating it
# in every sentinel-set panel would be redundant.
#
# Layout matches Fig. 1 (fig_ews_all_kc.jl): 3 columns (A=Complete,
# B=Scale-free, C=Bat), left to right, x 3 rows (row 1: R_S; row 2:
# n*Var_t R_S; row 3: microscopic). Only the K_c reference line is drawn
# (dotted); the 0.5*K_c reference used elsewhere is omitted here since S2/S3
# are already anchored to their own 0.9*K_c/0.5*K_c windows and marking a
# third point would clutter the panel. Legend (S2..S6 only) sits inside the
# first panel. Panel labels (A1)-(A3), (B1)-(B3), (C1)-(C3).
#
# Run (after simulate/01_ews_sentinel_sweep.jl):
#   julia --project=. figures/fig_ews_sentinel.jl
#
# Writes ../draft/ews_sentinel.pdf.

include(joinpath(@__DIR__, "..", "core", "kuramoto_core.jl"))   # defines SimParams, needed to deserialize bundle_ews_Kc20.jls
include(joinpath(@__DIR__, "..", "core", "figure_style.jl")); using .FigStyle
using CairoMakie, Serialization, Statistics, Printf
CairoMakie.update_theme!(Axis = (; xgridvisible = false, ygridvisible = false))

const HERE  = @__DIR__
const CACHE = joinpath(HERE, "..", "cache")
const OUTDIR = normpath(joinpath(HERE, "..", "draft")); mkpath(OUTDIR)
const B = deserialize(joinpath(CACHE, "bundle_ews_Kc20.jls"))
const KEYS = ["complete", "scale_free", "animal_18_1"]   # Complete, Scale-free, Bat (left to right, as in Fig. 1)
const REALIZATION = 1
const RUN = 2   # single representative run, shown raw (see header note)
const ROWS = [findfirst(m -> String(m.key) == k, B.meta) for k in KEYS]

# ---- style knobs ----
const NNET = length(ROWS)
const FIGSIZE = (2100, 1600)
const ROWLABELSIZE = 44
const PANELLABELSIZE = 40
const TICKLABELSIZE = 32
const SET_LABELS = ["S1", "S2", "S3", "S4", "S5", "S6"]   # legend needs only the short set names
# S1's slot (black) is kept reserved but unused, so S2..S6 keep the same
# colour/marker identity as in any other figure that does include S1:
const PAL = [RGBf(0, 0, 0), RGBf(0.0, 0.45, 0.70), RGBf(0.90, 0.45, 0.0), RGBf(0.0, 0.60, 0.50),
             RGBf(0.58, 0.20, 0.60), RGBf(0.85, 0.10, 0.10)]
const MARKERS = [:circle, :rect, :utriangle, :diamond, :star5, :xcross]
const OTHER_LW = 4.
const MARKERSIZE = 30; const N_MARKS = 10
const KC_LW = 3.0

"""Bundle column (1..8, S1/S2-1/S2-2/S3-1/S3-2/S4/S5/S6 -- see
simulate/02's header note) for sentinel set `newset` in 2:6 and `signal` in
(:RS, :RSvar, :micro). S2/S3 resolve to the c_i-based column (S2-1/S3-1) for
R_S and Var_t R_S, but to the Var_t delta_theta-based column (S2-2/S3-2) for
the microscopic EWS; S4 (hub), S5 (periphery), and S6 (random) resolve to a
single column regardless of signal -- see main.tex Sec. 'Selection of
sentinel nodes'."""
function newset_column(newset, signal)
    newset == 2 && return signal == :micro ? 3 : 2
    newset == 3 && return signal == :micro ? 5 : 4
    newset == 4 && return 6   # S4: hub
    newset == 5 && return 7   # S5: periphery
    newset == 6 && return 8   # S6: random
    error("newset must be 2..6 (S1 is handled separately, see fig_ews_all_kc.jl)")
end

"Column c of `field` for run `rr` (from bundle `B`), for sentinel set
`newset` (2..6) and `signal`. The macroscopic field (RSvar) is rescaled by
the set size n to make it intensive (Var_t R_S ~ 1/n), per main.tex Sec.
'Early warning signals'."
function set_column(rr, field, newset, signal)
    col = newset_column(newset, signal)
    values = getfield(rr, field)[:, col]
    field === :RSvar ? values .* length(rr.sets[col]) : values
end

"`Krel` is K/K_c; K_c itself is always at 1 once rescaled. `rr` is the
single representative run's record (raw trajectory, not averaged across
runs -- see header note)."
function ews_panel!(ax, rr, Krel, field, signal; show_legend = false)
    for newset in 2:6                     # S1 (all nodes) omitted -- see fig_ews_all_kc.jl instead
        vals = set_column(rr, field, newset, signal)
        lines!(ax, Krel, vals; color = PAL[newset], linewidth = OTHER_LW, label = SET_LABELS[newset])
    end
    vlines!(ax, [1.0]; color = :black, linestyle = :dot, linewidth = KC_LW)
    xlims!(ax, 0, 2.0)
    show_legend && axislegend(ax; position = :rb, framevisible = false, labelsize = FigStyle.LEGENDSIZE-5,
                               nbanks = 1, merge = true, patchsize = (70, 3))
end

row_letter(row) = string(Char('A' + row - 1))

# (field, signal, ylabel)
const COLS_SPEC = [
    (:RS,    :RS,    L"\langle R_S \rangle_t"),
    (:RSvar, :RSvar, L"n \,\, \mathrm{Var}(R_S)"),
    (:SM,    :micro, L"V_S"),
]

fig = Figure(size = FIGSIZE, figure_padding = (20, 30, 10, 10))

for (row, i) in enumerate(ROWS)
    m = B.meta[i]; rm = B.realization_meta[i][REALIZATION]
    rr = B.results_allfreq[i][REALIZATION][RUN]
    Krel = rr.Kv ./ rm.Kc

    Label(fig[row, 0], m.label; rotation = pi/2, tellheight = false, fontsize = ROWLABELSIZE)

    for (col, (field, signal, ylab)) in enumerate(COLS_SPEC)
        ax = FigStyle.styled_axis(fig[row, col];
            xlabel = row == length(COLS_SPEC) ? L"K/K_\mathrm{c}" : "", ylabel = ylab,
            xticklabelsvisible = row == length(COLS_SPEC),
            xticklabelsize = TICKLABELSIZE, yticklabelsize = TICKLABELSIZE,
            yticks = col == 1 ? collect(0.25:0.25:1.0) : Makie.automatic)
        Label(fig[row, col, TopLeft()], "($(row_letter(row))$(col))";
              fontsize = PANELLABELSIZE, font = :bold, padding = (80, -80, 5, 0), halign = :left)
        ews_panel!(ax, rr, Krel, field, signal; show_legend = (col == 1 && row == 1))
        # All three rows: shared range across all 9 networks (this figure +
        # fig_ews_sentinel_appendix.jl, across S2-S6 curves).
        col == 1 && ylims!(ax, 0.3, 1.05)
        col == 2 && ylims!(ax, -0.01, 0.45)
        col == 3 && ylims!(ax, -0.5, 19.5)
    end
end

colsize!(fig.layout, 0, Fixed(40))
for col in 1:NNET
    colsize!(fig.layout, col, Relative(1 / NNET))
end
colgap!(fig.layout, 50)
rowgap!(fig.layout, 30)

save(joinpath(OUTDIR, "ews_sentinel_rev.pdf"), fig)
println("wrote ews_sentinel_rev.pdf to ", OUTDIR)
