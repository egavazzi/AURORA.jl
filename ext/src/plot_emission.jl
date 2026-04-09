using Makie
using MAT: matread

# ======================================================================================== #
# Single-axis functions (plot_emission!, plot_column_emission!)
# ======================================================================================== #

"""
    _height_of_max(Q, h_atm)

Compute the altitude of maximum emission/ionization for each time step.
Values are set to NaN where the maximum is less than 1/10 of the global maximum.
"""
function _height_of_max(Q, h_atm)
    h_max = [h_atm[i_max[1]] for i_max in vec(findmax(Q, dims=1)[2])]
    h_max[vec(maximum(Q, dims=1)) .< maximum(Q) / 10] .= NaN
    return h_max
end

const _EMISSION_LABELS = Dict(
    :Q4278 => "4278 Å",
    :Q6730 => "6730 Å",
    :Q7774 => "7774 Å",
    :Q8446 => "8446 Å",
    :QO1D  => "O(¹D)",
    :QO1S  => "O(¹S)",
)

function AURORA.plot_emission!(ax, data::AURORA.VolumeExcitationResult, wavelength::Symbol;
                               colorrange = nothing, colormap = :viridis,
                               show_height_of_max = true, kwargs...)
    Q = getfield(data, wavelength)
    h_km = data.h_atm ./ 1e3
    t = data.t

    if colorrange === nothing
        Qmax = maximum(Q)
        colorrange = (max(Qmax / 1e4, 1.0), Qmax)
    end

    hm = Makie.heatmap!(ax, t, h_km, Q'; colorscale = log10, colorrange, colormap,
                        rasterize = true, kwargs...)

    if show_height_of_max
        h_max = _height_of_max(Q, h_km)
        Makie.lines!(ax, t, h_max; color = :red, linestyle = :dash, linewidth = 2)
    end

    return hm
end

const _COLUMN_WAVELENGTHS = [
    (:I_4278, Makie.rich("I", Makie.subscript("4278")), :blue),
    (:I_6730, Makie.rich("I", Makie.subscript("6730")), :red),
    (:I_7774, Makie.rich("I", Makie.subscript("7774")), Makie.RGBf(0.5, 0, 0)),
    (:I_8446, Makie.rich("I", Makie.subscript("8446")), :black),
    (:I_O1D,  Makie.rich("q", Makie.subscript("O(¹D)")), Makie.RGBf(1, 0.2, 0)),
    (:I_O1S,  Makie.rich("q", Makie.subscript("O(¹S)")), :green),
]

function AURORA.plot_column_emission!(ax, data::AURORA.ColumnExcitationResult;
                                      wavelengths = [:I_4278, :I_6730, :I_7774, :I_8446, :I_O1D, :I_O1S],
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


# ======================================================================================== #
# Full-figure functions (plot_emission, plot_column_emission)
# ======================================================================================== #

function AURORA.plot_emission(data::AURORA.VolumeExcitationResult;
                              size = (1000, 800), colorrange = nothing, kwargs...)
    fig = Makie.Figure(; size, fontsize = 20)

    ax_kwargs = (; xticksmirrored = true, yticksmirrored = true,
                   xminorticksvisible = true, yminorticksvisible = true,
                   limits = ((0, data.t[end]), (100, 400)))

    wavelengths = [:Q4278, :Q6730, :Q7774, :Q8446]
    positions = [(1, 1), (1, 2), (2, 1), (2, 2)]

    for (wl, (r, c)) in zip(wavelengths, positions)
        g = fig[r, c] = Makie.GridLayout()
        extra_kwargs = (;)
        if r == 1
            extra_kwargs = (; xticklabelsvisible = false,)
        end
        if c == 2
            extra_kwargs = (; extra_kwargs..., yticklabelsvisible = false,)
        end
        xlabel = r == 2 ? "t (s)" : ""
        ylabel = c == 1 ? "Altitude (km)" : ""
        ax = Makie.Axis(g[1, 1]; ax_kwargs..., title = _EMISSION_LABELS[wl],
                        xlabel, ylabel, extra_kwargs...)
        hm = AURORA.plot_emission!(ax, data, wl; colorrange, kwargs...)
        Makie.Colorbar(g[1, 2], hm; label = "photons/m³/s")
        Makie.colgap!(g, 10)
    end

    return fig
end

function AURORA.plot_column_emission(data::AURORA.ColumnExcitationResult;
                                     size = (1000, 800), kwargs...)
    fig = Makie.Figure(; size, fontsize = 20)
    ax = Makie.Axis(fig[1, 1]; xlabel = "t (s)", ylabel = "Intensity (R)",
                    yscale = log10,
                    xminorticksvisible = true, xminorgridvisible = true,
                    xminorticks = Makie.IntervalsBetween(10),
                    yminorticksvisible = true, yminorgridvisible = true,
                    yminorticks = Makie.IntervalsBetween(9), yticksmirrored = true)
    AURORA.plot_column_emission!(ax, data; kwargs...)
    Makie.Legend(fig[1, 2], ax; patchsize = [40, 20])

    return fig
end
