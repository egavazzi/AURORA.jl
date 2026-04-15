using Makie

function AURORA.plot_excitation!(ax, data::AURORA.VolumeExcitationResult;
                                 field::Symbol = :total,
                                 time_index::Union{Int,Nothing} = nothing,
                                 colorrange = nothing, colormap = :viridis,
                                 show_contours = false, show_height_of_max = true,
                                 kwargs...)
    Q = if field === :total
        data.QN2i .+ data.QO2i .+ data.QOi
    else
        getfield(data, field)
    end

    h_km = data.h_atm ./ 1e3
    t = data.t
    n_t = length(t)

    if n_t == 1 || time_index !== nothing
        i_t, _ = _select_time_index(t, n_t, time_index)
        q = Q[:, i_t]
        p = lines!(ax, q, h_km; kwargs...)
        xmin, xmax = _profile_xlims_from_top(q)
        xlims!(ax, xmin, xmax)
        if ax.xscale[] === log10
            ax.xticks = LogTicks(ceil(Int, log10(xmin)):floor(Int, log10(xmax)))
            ax.xminorticks = IntervalsBetween(9)
            ax.xminorticksvisible = true
        end
        return p
    end

    if colorrange === nothing
        Qmax = maximum(Q)
        colorrange = (max(Qmax / 1e4, 1.0), Qmax)
    end

    hm = heatmap!(ax, t, h_km, Q'; colorscale = log10, colorrange, colormap,
                  rasterize = true, kwargs...)

    if show_contours
        Q_log = log10.(Q)
        replace!(Q_log, -Inf => 0.0)
        contour!(ax, t, h_km, Q_log'; levels = 6:10, color = :black, labels = true)
    end

    if show_height_of_max
        h_max = _height_of_max(Q, h_km)
        lines!(ax, t, h_max; color = :red, linestyle = :dash, linewidth = 2)
    end

    return hm
end

const _COLUMN_WAVELENGTHS = [
    (:I_4278, Makie.rich("I", Makie.subscript("4278")), :blue),
    (:I_6730, Makie.rich("I", Makie.subscript("6730")), :red),
    (:I_7774, Makie.rich("I", Makie.subscript("7774")), Makie.RGBf(0.5, 0, 0)),
    (:I_8446, Makie.rich("I", Makie.subscript("8446")), :black),
    # Here we label O1S and O1D as excitation/production rates (q) to avoid confusion with
    # the quasi-instantaneous emission rates (I) of the other wavelengths
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
    lines_plots = Makie.Lines[]

    for (field, label, color) in _COLUMN_WAVELENGTHS
        field in wavelengths || continue
        I = getfield(data, field) ./ 1e10  # convert to Rayleigh
        linestyle = field in (:I_O1D, :I_O1S) ? :dash : :solid
        p = Makie.lines!(ax, t, I; label, color, linestyle, kwargs...)
        push!(lines_plots, p)
        local_max = maximum(filter(isfinite, I); init = -Inf)
        ymax = max(ymax, local_max)
    end

    # Set ylims manually to work around Makie log10-scale autolimit issues
    if isfinite(ymax) && ymax > 0
        Makie.ylims!(ax, ymax / 1e3, ymax * 2)
    end

    return lines_plots
end
