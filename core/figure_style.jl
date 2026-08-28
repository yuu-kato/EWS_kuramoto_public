# figure_style.jl
#
# Shared CairoMakie style constants for the figures in figures/. Keeping these
# in one place means every figure uses the same font sizes and axis
# conventions; tweak the numbers here (or override per-figure) rather than
# hand-editing each figure script.
#
# Usage:
#     include("../core/figure_style.jl"); using .FigStyle
#     fig = Figure(size = (1000, 620))
#     ax  = FigStyle.styled_axis(fig[1, 1]; xlabel = "coupling K", ylabel = "R")

module FigStyle

using CairoMakie

const LABELSIZE  = 40   # axis label size
const TICKSIZE   = 30   # tick label size
const TITLESIZE  = 40   # panel title size (use ~28-34 for long/multi-word titles)
const LEGENDSIZE = 40   # legend label size

"Axis with the canonical label/tick/title sizes pre-applied. Pass
xticklabelsize/yticklabelsize to override the tick label size for a single
figure without changing the shared TICKSIZE default."
function styled_axis(pos; xlabel = "", ylabel = "", title = "",
                      titlesize = TITLESIZE,
                      xticklabelsize = TICKSIZE, yticklabelsize = TICKSIZE, kwargs...)
    Axis(pos; xlabel = xlabel, ylabel = ylabel, title = title,
         xlabelsize = LABELSIZE, ylabelsize = LABELSIZE,
         xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize,
         titlesize = titlesize, kwargs...)
end

export LABELSIZE, TICKSIZE, TITLESIZE, LEGENDSIZE, styled_axis

end # module
