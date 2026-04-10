using Makie

function AURORA.plot_input!(ax, data::AURORA.IeTopResult;
                            beams = 1,
                            colorrange = nothing,
                            colormap = :inferno,
                            kwargs...)
    all_weights = AURORA.beam_weight(acosd.(data.mu_lims))
    beams = beams isa Int ? [beams] : beams # always treat as array for indexing
    BeamW = sum(all_weights[beams])
    data_heatmap = dropdims(sum(data.Ietop[beams, :, :]; dims=1); dims=1) ./ BeamW ./ data.ΔE' .* data.E_centers'

    if colorrange === nothing
        dmax = maximum(data_heatmap)
        colorrange = (dmax / 1e4, dmax)
    end

    return heatmap!(ax, data.t, data.E_centers, data_heatmap;
                    colorscale = log10, colorrange, colormap,
                    rasterize = true, kwargs...)
end

function AURORA.plot_input(sim::AURORA.AuroraSimulation)
    model = sim.model
    E_centers = model.energy_grid.E_centers
    ΔE = model.energy_grid.ΔE
    θ_lims = model.pitch_angle_grid.θ_lims
    Ω_beam = model.scattering.Ω_beam
    flux = sim.flux

    if isnothing(sim.time)
        # Steady-state: compute flux [n_beams, 1, n_E]
        Ie_top = AURORA.compute_flux(flux, model)
        return _plot_input_steady_state(Ie_top, E_centers, ΔE, θ_lims, Ω_beam, flux)
    else
        # Time-dependent: compute flux [n_beams, n_t, n_E]
        Ie_top = AURORA.compute_flux(flux, model, sim.time.t)
        t = collect(sim.time.t)
        return _plot_input_time_dependent(Ie_top, t, E_centers, ΔE, θ_lims, Ω_beam, flux)
    end
end


function _plot_input_steady_state(Ie_top, E_centers, ΔE, θ_lims, Ω_beam, flux)
    active_beams = flux.beams
    n_panels = length(active_beams)

    fig = Figure(size = (600, 400 * n_panels))

    for (idx, i_μ) in enumerate(active_beams)
        # Differential number flux: #e⁻/m²/s/eV/ster
        Ie_differential = Ie_top[i_μ, 1, :] ./ ΔE ./ Ω_beam[i_μ]

        θ_lo = round(min(θ_lims[i_μ], θ_lims[i_μ+1]); digits=1)
        θ_hi = round(max(θ_lims[i_μ], θ_lims[i_μ+1]); digits=1)
        ax = Axis(fig[idx, 1];
                  xlabel = "Energy (eV)",
                  ylabel = "Flux (#e⁻/m²/s/eV/ster)",
                  title = "Beam $i_μ (θ = $(θ_hi) – $(θ_lo)°)",
                  xscale = log10,
                  yscale = log10,
                  xminorticks = IntervalsBetween(9),
                  xminorticksvisible = true)
        lines!(ax, E_centers, Ie_differential)
    end

    Label(fig[0, :], "Input flux — Steady state"; fontsize = 16, tellwidth = false)

    return fig
end


function _plot_input_time_dependent(Ie_top, t, E_centers, ΔE, θ_lims, Ω_beam, flux)
    active_beams = flux.beams
    n_panels = length(active_beams)

    fig = Figure(size = (800, 600 * n_panels))

    for (idx, i_μ) in enumerate(active_beams)
        # Differential energy flux: eV/m²/s/eV/ster  (multiply number flux by E)
        IeE_diff = zeros(length(t), length(E_centers))
        for iE in eachindex(E_centers)
            IeE_diff[:, iE] .= Ie_top[i_μ, :, iE] .* E_centers[iE] ./ ΔE[iE] ./ Ω_beam[i_μ]
        end

        θ_lo = round(min(θ_lims[i_μ], θ_lims[i_μ+1]); digits=1)
        θ_hi = round(max(θ_lims[i_μ], θ_lims[i_μ+1]); digits=1)

        ax = Axis(fig[idx, 1];
                  xlabel = "Time (s)",
                  ylabel = "Energy (eV)",
                  title = "Beam $i_μ (θ = $(θ_hi) – $(θ_lo)°)",
                  yscale = log10,
                  yminorticks = IntervalsBetween(9),
                  yminorticksvisible = true)

        max_val = maximum(IeE_diff)
        if max_val > 0
            cr = (max_val / 1e4, max_val)
        else
            cr = (1e-1, 1e0)
        end

        hm = heatmap!(ax, t, E_centers, IeE_diff;
                      colorscale = log10,
                      colormap = :inferno,
                      colorrange = cr)
        Colorbar(fig[idx, 2], hm; label = "IeE (eV/m²/s/eV/ster)")
    end

    Label(fig[0, :], "Input flux — Time dependent"; fontsize = 16, tellwidth = false)

    return fig
end
