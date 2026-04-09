using Makie
using MAT: matread

# ======================================================================================== #
# Single-axis function (plot_ionization!)
# ======================================================================================== #

function AURORA.plot_ionization!(ax, data::AURORA.VolumeExcitationResult;
                                 species::Symbol = :total,
                                 colorrange = nothing, colormap = :viridis,
                                 show_contours = true, show_height_of_max = true,
                                 kwargs...)
    Q = if species === :total
        data.QN2i .+ data.QO2i .+ data.QOi
    else
        getfield(data, species)
    end

    h_km = data.h_atm ./ 1e3
    t = data.t

    if colorrange === nothing
        Qmax = maximum(Q)
        colorrange = (max(Qmax / 1e4, 1.0), Qmax)
    end

    hm = Makie.heatmap!(ax, t, h_km, Q'; colorscale = log10, colorrange, colormap,
                        rasterize = true, kwargs...)

    if show_contours
        Q_log = log10.(Q)
        replace!(Q_log, -Inf => 0.0)
        Makie.contour!(ax, t, h_km, Q_log'; levels = 6:10, color = :black, labels = true)
    end

    if show_height_of_max
        h_max = _height_of_max(Q, h_km)
        Makie.lines!(ax, t, h_max; color = :red, linestyle = :dash, linewidth = 2)
    end

    return hm
end


# ======================================================================================== #
# Full-figure function (plot_ionization)
# ======================================================================================== #

function AURORA.plot_ionization(data::AURORA.VolumeExcitationResult;
                                plot_Ietop::Bool = false,
                                angle_cone = [170, 180],
                                size = (850, 950), kwargs...)
    if plot_Ietop && data.savedir === nothing
        error("plot_Ietop=true requires data to have a known savedir. " *
              "Load data with load_volume_excitation(directory) to populate it.")
    end

    fig_height = plot_Ietop ? size[2] : 500
    fig = Makie.Figure(; size = (size[1], fig_height), fontsize = 20)

    row = 1
    ax_Ietop = nothing

    if plot_Ietop
        has_Ietop = any(startswith(f, "Ie_incoming_") for f in readdir(data.savedir))
        if !has_Ietop
            @warn "No Ie_incoming file found in $(data.savedir), skipping top-flux panel."
            plot_Ietop = false
        end
    end

    if plot_Ietop
        Ietop_file = AURORA.find_Ietop_file(data.savedir)
        ie_data = matread(Ietop_file)
        Ietop = ie_data["Ie_total"]
        μ_lims = ie_data["mu_lims"]
        t_top = ie_data["t_top"]
        t_top = [t_top; t_top[end] + diff(t_top)[end]] .- t_top[1]
        Egrid = vec(ie_data["E_centers"])
        dEgrid = vec(ie_data["dE"])

        ax_Ietop = Makie.Axis(fig[row, 1]; yscale = log10, ylabel = "Energy (eV)",
                              yminorticksvisible = true, yminorticks = Makie.IntervalsBetween(9),
                              xticklabelsvisible = false, xminorticksvisible = true,
                              xticksmirrored = true, yticksmirrored = true)
        θ_lims = acosd.(μ_lims)
        idx_θ = vec(angle_cone[1] .<= abs.(acosd.(AURORA.mu_avg(θ_lims))) .<= angle_cone[2])
        BeamW = AURORA.beam_weight([angle_cone[1], angle_cone[2]])
        data_heatmap = dropdims(sum(Ietop[idx_θ, :, :]; dims=1); dims=1) ./ BeamW ./ dEgrid' .* Egrid'
        hm_Ietop = Makie.heatmap!(ax_Ietop, t_top, Egrid, data_heatmap;
                                  colorrange = (1e6, maximum(data_heatmap)),
                                  colorscale = log10, colormap = :inferno)
        Makie.Colorbar(fig[row, 2], hm_Ietop; label = "IeE (eV/m²/s/eV/ster)")
        row += 1
    end

    h_km = data.h_atm ./ 1e3
    t = data.t

    ax_Q = Makie.Axis(fig[row, 1]; xlabel = "t (s)", ylabel = "Altitude (km)",
                      xminorticksvisible = true, yminorticksvisible = true,
                      xticksmirrored = true, yticksmirrored = true)
    hm_Q = AURORA.plot_ionization!(ax_Q, data; kwargs...)
    Makie.Colorbar(fig[row, 2], hm_Q; label = "Ionization rate (/m³/s)")

    Makie.colsize!(fig.layout, 1, Makie.Auto())
    if plot_Ietop
        Makie.rowsize!(fig.layout, 1, Makie.Relative(1/5))
        Makie.linkxaxes!(ax_Ietop, ax_Q)
        Makie.xlims!(ax_Q, 0, min(t_top[end], t[end]))
    end
    Makie.ylims!(ax_Q, 100, 500)

    return fig
end
