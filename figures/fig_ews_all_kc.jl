# figures/fig_ews_all_kc.jl
#
# Fig. 1 (main.tex, label fig:kc): critical-coupling determination for the
# three representative networks shown in the main text -- Complete,
# Scale-free, and Bat (the other six networks are shown in
# Appendix figures/fig_ews_all_kc_appendix.jl). Illustrative figure: uses only
# the first of the NFREQ frequency realizations' second run
# (cache/bundle_ews_Kc20.jls's `results_allfreq[i][1][2]`).
#
# This figure's data comes entirely from simulate/01_ews_sentinel_sweep.jl's
# output bundle: that bundle's S1 sentinel set ("all N nodes", bundle column
# 1) IS the full-observation quantities this figure needs -- <R>_t, Var_t R,
# and the microscopic EWS -- computed over exactly the K-range (0.01 to
# 2*Kc20) this figure plots (K/K_c in [0,2]), so no separate simulation is
# needed here. See simulate/01_ews_sentinel_sweep.jl's header note for the
# K_c-determination-first, EWS-calculation-second rationale.
#
# Layout: 3 columns (A=Complete, B=Scale-free, C=Bat) x 3 rows (row 1: global
# order parameter <R>_t; row 2: rescaled macroscopic variance N*Var_t R(t);
# row 3: microscopic EWS mean_i Var_t[delta_theta_i] under full observation
# (S1)); both rows 2 and 3 are UNBIASED temporal variances --
# core/kuramoto_core.jl's `var_unbiased`/`resid_var`. All three rows show a
# single representative run (run 2), raw and unaveraged, of the first
# frequency realization. X-axis is K/K_c (each network's own, realization-1
# K_c), 0 to 2. Dotted black: K_c (always at 1; core/kuramoto_core.jl's
# `Kc_threshold`, main.tex Eq. def_kcp), the sole anchor used throughout the
# paper. Panel labels (A1)-(A3), (B1)-(B3), (C1)-(C3).
#
# Run (after simulate/01_ews_sentinel_sweep.jl):
#   julia --project=. figures/fig_ews_all_kc.jl
#
# Writes ../draft/ews_all_kc.pdf.

include(joinpath(@__DIR__, "..", "core", "kuramoto_core.jl"))   # defines SimParams, needed to deserialize bundle_ews_Kc20.jls
include(joinpath(@__DIR__, "..", "core", "figure_style.jl")); using .FigStyle
using CairoMakie, Serialization, Statistics, Printf
CairoMakie.update_theme!(Axis = (; xgridvisible = false, ygridvisible = false))

const HERE  = @__DIR__
const CACHE = joinpath(HERE, "..", "cache")
const OUTDIR = normpath(joinpath(HERE, "..", "draft")); mkpath(OUTDIR)
const B = deserialize(joinpath(CACHE, "bundle_ews_Kc20.jls"))
const S1_COL = 1   # bundle column for S1 ("all N nodes") -- see simulate/01_ews_sentinel_sweep.jl's header note
const KEYS = ["complete", "scale_free", "animal_18_1"]   # Complete, Scale-free, Bat
const REALIZATION = 1
const RUN = 2   # single representative run (see header note)
const ROWS = [findfirst(m -> String(m.key) == k, B.meta) for k in KEYS]

const FIGSIZE = (2100, 1600)
const ROWLABELSIZE = 44
const PANELLABELSIZE = 40
const TICKLABELSIZE = 32
const SINGLE_LW = 5.0
const KC_LW = 3.0
const KC20_COLOR = :black

row_letter(row) = string(Char('A' + row - 1))

fig = Figure(size = FIGSIZE, figure_padding = (20, 30, 10, 10))

for (row, i) in enumerate(ROWS)
    m = B.meta[i]; rm = B.realization_meta[i][REALIZATION]
    rr = B.results_allfreq[i][REALIZATION][RUN]
    Krel = rr.Kv ./ rm.Kc
    N = m.N

    Label(fig[row, 0], m.label; rotation = pi/2, tellheight = false, fontsize = ROWLABELSIZE)

    # ---- row 1: global order parameter <R>_t, single run ----
    ax1 = FigStyle.styled_axis(fig[row, 1]; xlabel = row == 3 ? L"K/K_\mathrm{c}" : "", ylabel = L"\langle R \rangle_t",
        yticklabelsize = TICKLABELSIZE, xticklabelsvisible = row == 3, xticks = 0:0.5:2.0, yticks = 0:0.25:1.0)
    Label(fig[row, 1, TopLeft()], "($(row_letter(row))1)"; fontsize = PANELLABELSIZE, font = :bold, padding = (110, -110, 15, 0), halign = :left)
    lines!(ax1, Krel, rr.RS[:, S1_COL]; color = :navy, linewidth = SINGLE_LW)
    vlines!(ax1, [1.0]; color = KC20_COLOR, linestyle = :dot, linewidth = KC_LW, label = L"K_\mathrm{c}")
    row == 1 && axislegend(ax1; position = :lt, framevisible = false, labelsize = FigStyle.LEGENDSIZE, patchsize = (60, 3))
    xlims!(ax1, 0, 2.0)
    ylims!(ax1, 0.0, 1.0)

    # ---- row 2: rescaled macroscopic variance N.Var_t R(t), single run ----
    # y-range shared across all 9 networks (this figure + fig_ews_all_kc_appendix.jl).
    ax2 = FigStyle.styled_axis(fig[row, 2]; xlabel = row == 3 ? L"K/K_\mathrm{c}" : "", ylabel = L"N \,\, \mathrm{Var}(R)",
        yticklabelsize = TICKLABELSIZE,
        xticklabelsvisible = row == 3, xticks = 0:0.5:2.0)
    Label(fig[row, 2, TopLeft()], "($(row_letter(row))2)"; fontsize = PANELLABELSIZE, font = :bold, padding = (80, -80, 15, 0), halign = :left)
    lines!(ax2, Krel, N .* rr.RSvar[:, S1_COL]; color = :firebrick, linewidth = SINGLE_LW)
    vlines!(ax2, [1.0]; color = KC20_COLOR, linestyle = :dot, linewidth = KC_LW)
    xlims!(ax2, 0, 2.0)
    ylims!(ax2, -0.05, 1.7)

    # ---- row 3: microscopic EWS, full observation (S1), single run ----
    # y-range shared across all 9 networks (this figure + appendix); ticks fixed at 0,2,4,6.
    ax3 = FigStyle.styled_axis(fig[row, 3]; xlabel = row == 3 ? L"K/K_\mathrm{c}" : "",
        ylabel = L"V_S",
        xticklabelsize = TICKLABELSIZE, xticklabelsvisible = row == 3, yticklabelsize = TICKLABELSIZE, xticks = 0:0.5:2.0, yticks = 0:2:8)
    Label(fig[row, 3, TopLeft()], "($(row_letter(row))3)"; fontsize = PANELLABELSIZE, font = :bold, padding = (60, -60, 15, 0), halign = :left)
    lines!(ax3, Krel, rr.SM[:, S1_COL]; color = :darkgreen, linewidth = SINGLE_LW)
    vlines!(ax3, [1.0]; color = KC20_COLOR, linestyle = :dot, linewidth = KC_LW)
    xlims!(ax3, 0, 2.0)
    ylims!(ax3, -0.1, 9.0)
end

colsize!(fig.layout, 0, Fixed(40))
for col in 1:3
    colsize!(fig.layout, col, Relative(1 / length(KEYS)))
end
colgap!(fig.layout, 50)
rowgap!(fig.layout, 30)

save(joinpath(OUTDIR, "ews_all_kc.pdf"), fig)
println("wrote ews_all_kc.pdf to ", OUTDIR)
