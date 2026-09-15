using Makie
using Printf

# Where the input energy flux ends up, bottom to top of the stack: field name, legend name,
# colour, and the colour its in-bar label needs to stay readable.
const _ENERGY_BUDGET_SEGMENTS = (
    (:ionization, "ionization", Makie.RGBf(0.16, 0.44, 0.30), :white),
    (:excitation, "excitation (non-ionizing)", Makie.RGBf(0.42, 0.72, 0.50), :white),
    (:heating, "thermal e-e loss", Makie.RGBf(0.45, 0.42, 0.72), :white),
    (:escape, "backscattered", Makie.RGBf(0.90, 0.60, 0.0), :black),
)

const _ENERGY_BUDGET_SEGMENT_COLORS = [c for (_, _, c, _) in _ENERGY_BUDGET_SEGMENTS]
const _ENERGY_BUDGET_SEGMENT_NAMES = [n for (_, n, _, _) in _ENERGY_BUDGET_SEGMENTS]

# Minimum spacing between two segment labels, as a fraction of the axis height. A segment
# thinner than this has its label drawn over the neighbouring segments.
const _SEGMENT_LABEL_MIN_GAP = 0.04

# Power of ten that puts the bar height between 1 and 10.
_budget_scale(input) = input > 0 ? 10.0^floor(log10(input)) : 1.0

const _SUPERSCRIPTS = Dict('0' => '⁰', '1' => '¹', '2' => '²', '3' => '³', '4' => '⁴',
                           '5' => '⁵', '6' => '⁶', '7' => '⁷', '8' => '⁸', '9' => '⁹',
                           '-' => '⁻')
_superscript(n::Integer) = String([_SUPERSCRIPTS[c] for c in string(n)])

# Scaled values, with the same guarantee: a nonzero term never reads as "0.0".
function _value_string(v)
    v != 0 && abs(v) < 0.005 && return string(round(v; sigdigits = 2))
    return string(round(v; digits = 2))
end

# Keep the axis ticks of every bar already drawn, so repeated `plot_energy_budget!` calls on
# one axis accumulate labelled positions instead of overwriting each other.
function _add_xtick!(ax, x, label)
    ticks = ax.xticks[]
    positions, labels = if ticks isa Tuple
        Float64.(collect(ticks[1])), String.(collect(ticks[2]))
    else
        Float64[], String[]
    end
    i = findfirst(==(Float64(x)), positions)
    if i === nothing
        push!(positions, Float64(x))
        push!(labels, label)
    else
        labels[i] = label
    end
    order = sortperm(positions)
    ax.xticks = (positions[order], labels[order])
    return nothing
end

function AURORA.plot_energy_budget!(ax, budget::AURORA.EnergyBudget;
                                    x = 1, label = nothing, scale = nothing)
    input = budget.input
    sc = something(scale, _budget_scale(input))
    values = [getproperty(budget, field) / sc for (field, _, _, _) in _ENERGY_BUDGET_SEGMENTS]
    top = max(input / sc, sum(values))

    plot = barplot!(ax, fill(Float64(x), length(values)), values;
                    stack = eachindex(values), color = _ENERGY_BUDGET_SEGMENT_COLORS)
    hlines!(ax, [input / sc]; color = :black, linestyle = :dash, linewidth = 2)

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
        text = "$(_value_string(values[i])) " *
               "($(AURORA.percent_string(getproperty(budget, field), input))%)"
        text!(ax, Float64(x), heights[i]; text, align = (:center, :center), fontsize = 14,
              color = textcolor)
    end

    label !== nothing && _add_xtick!(ax, x, label)
    return plot
end

function AURORA.plot_energy_budget(budget::AURORA.EnergyBudget; label = nothing)
    input = budget.input
    sc = _budget_scale(input)
    exponent = round(Int, log10(sc))
    accounted = budget.inelastic + budget.heating + budget.escape
    ratio = input != 0 ? accounted / input : NaN

    fig = Figure(size = (560, 760), fontsize = 17)
    ax = Axis(fig[1, 1];
              ylabel = "Vertical energy flux (×10$(_superscript(exponent)) " *
                       "$(AURORA.energy_units(budget)))",
              title = @sprintf("(deposited + backscattered) / input = %.3f", ratio),
              titlefont = :regular, titlesize = 14,
              yminorticksvisible = true,
              xticks = ([1.0], [something(label, "run")]))
    AURORA.plot_energy_budget!(ax, budget; x = 1, scale = sc)

    xmin, xmax = 0.3, 1.7
    text!(ax, xmin + 0.01 * (xmax - xmin), input / sc; text = "input",
          align = (:left, :bottom), fontsize = 13, color = :black)
    Legend(fig[2, 1],
           [PolyElement(color = c) for c in _ENERGY_BUDGET_SEGMENT_COLORS],
           _ENERGY_BUDGET_SEGMENT_NAMES;
           orientation = :horizontal, framevisible = true, nbanks = 2)

    stack = sum(getproperty(budget, field) for (field, _, _, _) in _ENERGY_BUDGET_SEGMENTS)
    ylims!(ax, 0, 1.15 * max(input, stack) / sc)
    xlims!(ax, xmin, xmax)
    return fig
end
