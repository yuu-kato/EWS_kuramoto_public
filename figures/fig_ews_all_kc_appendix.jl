# figures/fig_ews_all_kc_appendix.jl
#
# Appendix counterpart to figures/fig_ews_all_kc.jl (main.tex Fig. 1): the
# same critical-coupling-determination panels, but for the six networks NOT
# shown in the main text -- ER, Highschool boy, House finch, Surfer, Hall,
# Jazz. Illustrative figure: uses only the first of the NFREQ frequency
# realizations' second run (cache/bundle_ews_Kc20.jls's
# `results_allfreq[i][1][2]`) -- see figures/fig_ews_all_kc.jl's header note
# for why this figure's data comes entirely from
# simulate/01_ews_sentinel_sweep.jl's bundle (S1 = "all N nodes"), rather
# than a dedicated sweep of its own. Both the order-parameter variance
# (col 2) and the microscopic EWS (col 3) are
# UNBIASED temporal variances -- see core/kuramoto_core.jl's
# `var_unbiased`/`resid_var`.
#
# Layout: 6 rows (one per network) x 3 columns (col 1: global order
# parameter <R>_t; col 2: rescaled macroscopic variance N*Var_t R(t); col 3:
# microscopic EWS mean_i Var_t[delta_theta_i] under full observation (S1)).
# All three columns show a single representative run (run 2), raw and
# unaveraged, of the first frequency realization. X-axis is K/K_c (each
# network's own, realization-1 K_c), 0 to
# 2. Dashed orange: K_c (always at 1; core/kuramoto_core.jl's
# `Kc_threshold`, main.tex Eq. def_kcp), the sole anchor used throughout the
# paper. Panel labels (A1)-(A3) .. (F1)-(F3).
#
# Run (after simulate/01_ews_sentinel_sweep.jl):
#   julia --project=. figures/fig_ews_all_kc_appendix.jl
#
# Writes ../draft/ews_all_kc_appendix.pdf.

include(joinpath(@__DIR__, "..", "core", "kuramoto_core.jl"))   # defines SimParams, needed to deserialize bundle_ews_Kc20.jls
include(joinpath(@__DIR__, "..", "core", "figure_style.jl")); using .FigStyle
using CairoMakie, Serialization, Statistics, Printf
CairoMakie.update_theme!(Axis = (; xgridvisible = false, ygridvisible = false))

const HERE  = @__DIR__
const CACHE = joinpath(HERE, "..", "cache")
const OUTDIR = normpath(joinpath(HERE, "..", "draft")); mkpath(OUTDIR)
const B = deserialize(joinpath(CACHE, "bundle_ews_Kc20.jls"))
const S1_COL = 1   # bundle column for S1 ("all N nodes") -- see simulate/01_ews_sentinel_sweep.jl's header note
const KEYS = ["erdos_renyi", "highschool_boys", "animal_6", "surfersb", "hall", "jazz"]
const REALIZATION = 1
const RUN = 2   # single representative run (see header note)
const ROWS = [findfirst(m -> String(m.key) == k, B.meta) for k in KEYS]

const NNET = length(ROWS)
const ROWHEIGHT = 420
const FIGSIZE = (2200, 80 + ROWHEIGHT * NNET)
const ROWLABELSIZE = 44
const PANELLABELSIZE = 40
const TICKLABELSIZE = 32
const SINGLE_LW = 4.0
const KC_LW = 3.0
const KC20_COLOR = :black

row_letter(row) = string(Char('A' + row - 1))

fig = Figure(size = FIGSIZE, figure_padding = (20, 20, 10, 10))

for (row, i) in enumerate(ROWS)
    m = B.meta[i]; rm = B.realization_meta[i][REALIZATION]
    rr = B.results_allfreq[i][REALIZATION][RUN]
    Krel = rr.Kv ./ rm.Kc
    N = m.N

    Label(fig[row, 0], m.label; rotation = pi/2, tellheight = false, fontsize = ROWLABELSIZE)

    ax1 = FigStyle.styled_axis(fig[row, 1];
        xlabel = row == NNET ? L"K/K_\mathrm{c}" : "", ylabel = L"\langle R \rangle_t",
        yticklabelsize = TICKLABELSIZE,
        xticklabelsize = TICKLABELSIZE, xticks = 0:0.5:2.0, yticks = 0:0.25:1.0,
        xticklabelsvisible = row == length(ROWS))
    Label(fig[row, 1, TopLeft()], "($(row_letter(row))1)"; fontsize = PANELLABELSIZE, font = :bold, padding = (0, 14, 15, 0), halign = :left)
    lines!(ax1, Krel, rr.RS[:, S1_COL]; color = :navy, linewidth = SINGLE_LW)
    vlines!(ax1, [1.0]; color = KC20_COLOR, linestyle = :dot, linewidth = KC_LW, label = L"K_\mathrm{c}")
    row == 1 && axislegend(ax1; position = :lt, framevisible = false, labelsize = FigStyle.LEGENDSIZE,
        patchsize = (50, 3))
    xlims!(ax1, 0, 2.0)
    ylims!(ax1, 0.0, 1.0)

    # y-range shared with fig_ews_all_kc.jl -- see that file's row-2 comment
    ax2 = FigStyle.styled_axis(fig[row, 2];
        xlabel = row == NNET ? L"K/K_\mathrm{c}" : "", ylabel = L"N \,\, \mathrm{Var}(R)",
        xticklabelsize = TICKLABELSIZE, yticklabelsize = TICKLABELSIZE, xticks = 0:0.5:2.0,
        xticklabelsvisible = row == length(ROWS))
    Label(fig[row, 2, TopLeft()], "($(row_letter(row))2)"; fontsize = PANELLABELSIZE, font = :bold, padding = (0, 14, 15, 0), halign = :left)
    lines!(ax2, Krel, N .* rr.RSvar[:, S1_COL]; color = :firebrick, linewidth = SINGLE_LW)
    vlines!(ax2, [1.0]; color = KC20_COLOR, linestyle = :dot, linewidth = KC_LW)
    xlims!(ax2, 0, 2.0)
    ylims!(ax2, -0.05, 2.5)

    # y-range shared with fig_ews_all_kc.jl -- see that file's row-3 comment
    ax3 = FigStyle.styled_axis(fig[row, 3];
        xlabel = row == NNET ? L"K/K_\mathrm{c}" : "", ylabel = L"V_S",
        xticklabelsize = TICKLABELSIZE, yticklabelsize = TICKLABELSIZE, xticks = 0:0.5:2.0, yticks = 0:1:3,
        xticklabelsvisible = row == length(ROWS))
    Label(fig[row, 3, TopLeft()], "($(row_letter(row))3)"; fontsize = PANELLABELSIZE, font = :bold, padding = (0, 14, 15, 0), halign = :left)
    lines!(ax3, Krel, rr.SM[:, S1_COL]; color = :darkgreen, linewidth = SINGLE_LW)
    vlines!(ax3, [1.0]; color = KC20_COLOR, linestyle = :dot, linewidth = KC_LW)
    xlims!(ax3, 0, 2.0)
    ylims!(ax3, -0.1, 2.5)
end

colsize!(fig.layout, 0, Fixed(30))
for col in 1:3
    colsize!(fig.layout, col, Relative(0.97 / 3))
end
colgap!(fig.layout, 50)
rowgap!(fig.layout, 20)

save(joinpath(OUTDIR, "ews_all_kc_appendix.pdf"), fig)
println("wrote ews_all_kc_appendix.pdf to ", OUTDIR)
