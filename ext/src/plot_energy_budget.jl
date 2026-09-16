using Makie
using Printf

# Where the input energy flux ends up, bottom to top of the stack: field name, legend name,
# colour, and the colour its label needs to stay readable on that colour.
const _ENERGY_BUDGET_SEGMENTS = (
    (:ionization, "ionization", Makie.RGBf(0.16, 0.44, 0.30), :white),
    (:excitation, "excitation (non-ionizing)", Makie.RGBf(0.42, 0.72, 0.50), :white),
    (:heating, "thermal e-e loss", Makie.RGBf(0.45, 0.42, 0.72), :white),
    (:bottom_escape, "escape ↓ bottom", Makie.RGBf(0.40, 0.31, 0.22), :white),
    (:escape, "backscattered ↑ top", Makie.RGBf(0.90, 0.60, 0.0), :black),
)

const _ENERGY_BUDGET_SEGMENT_COLORS = [c for (_, _, c, _) in _ENERGY_BUDGET_SEGMENTS]
const _ENERGY_BUDGET_SEGMENT_NAMES = [n for (_, n, _, _) in _ENERGY_BUDGET_SEGMENTS]

# Minimum spacing between two segment labels, as a fraction of the bar height. A segment
# thinner than this has its label drawn over the neighbouring segments.
const _SEGMENT_LABEL_MIN_GAP = 0.04

# A segment below this share of the input carries no information worth a label.
const _SEGMENT_LABEL_MIN_SHARE = 1e-6

# Half-width of a bar, and of the dashed input line drawn across it.
const _BAR_HALF_WIDTH = 0.4

# Gutter between the outermost bars and the edge of the axis. Kept small so the bars stay
# wide enough to hold a segment label, which is written over them in a contrasting colour
# and would fall on the figure background if it overflowed.
const _BAR_GUTTER = 0.15

# Three significant digits in scientific notation, e.g. 2.95×10¹⁶, as rich text.
function _sci_text(x)
    x == 0 && return Makie.rich("0")
    mantissa, exponent = split(@sprintf("%.2e", x), 'e')
    return Makie.rich(string(mantissa, "×10"), Makie.superscript(string(parse(Int, exponent))))
end

_budget_ylabel(budget) = budget.interval === nothing ?
    "Vertical energy flux (eV m⁻² s⁻¹)" : "Vertical energy (eV m⁻²)"

_budget_ratio(budget) =
    (budget.inelastic + budget.heating + budget.escape + budget.bottom_escape) / budget.input

function AURORA.plot_energy_budget!(ax, budget::AURORA.EnergyBudget; x = 1)
    values = [getproperty(budget, field) for (field, _, _, _) in _ENERGY_BUDGET_SEGMENTS]
    top = max(budget.input, sum(values))

    plot = barplot!(ax, fill(Float64(x), length(values)), values;
                    stack = eachindex(values), color = _ENERGY_BUDGET_SEGMENT_COLORS,
                    width = 2 * _BAR_HALF_WIDTH)
    linesegments!(ax, [Makie.Point2f(x - _BAR_HALF_WIDTH, budget.input),
                       Makie.Point2f(x + _BAR_HALF_WIDTH, budget.input)];
                  color = :black, linestyle = :dash, linewidth = 2)

    # Label heights, bottom to top, pushed apart where two segments meet too closely for
    # their labels to be read side by side.
    heights = Float64[]
    cumulative = 0.0
    for v in values
        y = cumulative + v / 2
        !isempty(heights) && (y = max(y, last(heights) + _SEGMENT_LABEL_MIN_GAP * top))
        push!(heights, y)
        cumulative += v
    end

    for (i, (field, _, _, textcolor)) in enumerate(_ENERGY_BUDGET_SEGMENTS)
        v = values[i]
        v < _SEGMENT_LABEL_MIN_SHARE * budget.input && continue
        text = Makie.rich(_sci_text(v), " ($(AURORA.percent_string(v, budget.input))%)")
        text!(ax, Float64(x), heights[i]; text, align = (:center, :center), fontsize = 14,
              color = textcolor)
    end
    return plot
end

function AURORA.plot_energy_budget(budget::AURORA.EnergyBudget; label = nothing)
    fig = Figure(size = (700, 760), fontsize = 17)
    ax = Axis(fig[1, 1];
              ylabel = _budget_ylabel(budget),
              title = @sprintf("(deposited + escaped) / input = %.3f",
                                      _budget_ratio(budget)),
              titlefont = :regular, titlesize = 14,
              yminorticksvisible = true,
              xticks = ([1.0], [something(label, "run")]))
    AURORA.plot_energy_budget!(ax, budget; x = 1)

    xmin, xmax = 1 - _BAR_HALF_WIDTH - _BAR_GUTTER, 1 + _BAR_HALF_WIDTH + _BAR_GUTTER
    text!(ax, 1 - _BAR_HALF_WIDTH, budget.input; text = "input",
          align = (:left, :bottom), fontsize = 13, color = :black)
    _budget_legend!(fig)
    _budget_limits!(ax, (budget,), xmin, xmax)
    return fig
end

function AURORA.plot_energy_budget(budgets::AbstractVector{<:AURORA.EnergyBudget};
                                   labels = nothing)
    isempty(budgets) && throw(ArgumentError("no budgets to plot"))
    labels === nothing || length(labels) == length(budgets) ||
        throw(DimensionMismatch("got $(length(labels)) labels for $(length(budgets)) " *
                                "budgets"))
    n = length(budgets)
    fig = Figure(size = (300 + 260 * n, 700), fontsize = 17)
    ax = Axis(fig[1, 1];
              ylabel = _budget_ylabel(first(budgets)),
              yminorticksvisible = true,
              xticks = labels === nothing ? (1.0:n, fill("", n)) :
                       (1.0:n, collect(String.(labels))))
    for (i, budget) in enumerate(budgets)
        AURORA.plot_energy_budget!(ax, budget; x = i)
    end
    _budget_legend!(fig)
    _budget_limits!(ax, budgets, 1 - _BAR_HALF_WIDTH - _BAR_GUTTER,
                    n + _BAR_HALF_WIDTH + _BAR_GUTTER)
    return fig
end

_budget_legend!(fig) =
    Legend(fig[2, 1],
           [PolyElement(color = c) for c in _ENERGY_BUDGET_SEGMENT_COLORS],
           _ENERGY_BUDGET_SEGMENT_NAMES;
           orientation = :horizontal, framevisible = true, nbanks = 2)

function _budget_limits!(ax, budgets, xmin, xmax)
    top = maximum(max(b.input, b.inelastic + b.heating + b.escape + b.bottom_escape)
                  for b in budgets)
    ylims!(ax, 0, 1.15 * top)
    xlims!(ax, xmin, xmax)
    return nothing
end
