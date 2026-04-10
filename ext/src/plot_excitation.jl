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
        ax.xticks = LogTicks(ceil(Int, log10(xmin)):floor(Int, log10(xmax)))
        ax.xminorticks = IntervalsBetween(9)
        ax.xminorticksvisible = true
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
