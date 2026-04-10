using Makie

const _COLUMN_WAVELENGTHS = [
    (:I_4278, Makie.rich("I", Makie.subscript("4278")), :blue),
    (:I_6730, Makie.rich("I", Makie.subscript("6730")), :red),
    (:I_7774, Makie.rich("I", Makie.subscript("7774")), Makie.RGBf(0.5, 0, 0)),
    (:I_8446, Makie.rich("I", Makie.subscript("8446")), :black),
    (:I_O1D,  Makie.rich("q", Makie.subscript("O(¹D)")), Makie.RGBf(1, 0.2, 0)),
    (:I_O1S,  Makie.rich("q", Makie.subscript("O(¹S)")), :green),
]

function AURORA.plot_column_excitation!(ax, data::AURORA.ColumnExcitationResult;
                                        wavelengths = [:I_4278,
                                            :I_6730,
                                            :I_7774,
                                            :I_8446,
                                            :I_O1D,
                                            :I_O1S],
                                        kwargs...)
    t = data.t
    ymax = -Inf

    for (field, label, color) in _COLUMN_WAVELENGTHS
        field in wavelengths || continue
        I = getfield(data, field) ./ 1e10  # convert to Rayleigh
        linestyle = field in (:I_O1D, :I_O1S) ? :dash : :solid
        Makie.lines!(ax, t, I; label, color, linestyle, kwargs...)
        local_max = maximum(filter(isfinite, I); init = -Inf)
        ymax = max(ymax, local_max)
    end

    # Set ylims manually to work around Makie log10-scale autolimit issues
    if isfinite(ymax) && ymax > 0
        Makie.ylims!(ax, ymax / 1e3, ymax * 2)
    end

    return ax
end
