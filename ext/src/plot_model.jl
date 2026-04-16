using Makie

const ALL_PANELS = [:atmosphere, :energy_levels, :energy_grid, :cross_sections,
                    :phase_functions, :scattering]

function AURORA.plot_model(model::AURORA.AuroraModel; panels=[:all])
    if :all in panels
        panels = ALL_PANELS
    end

    figs = Dict{Symbol, Makie.Figure}()

    for panel in panels
        if panel === :atmosphere
            figs[:atmosphere] = _plot_atmosphere(model)
        elseif panel === :energy_levels
            figs[:energy_levels] = _plot_energy_levels(model)
        elseif panel === :energy_grid
            figs[:energy_grid] = _plot_energy_grid(model)
        elseif panel === :cross_sections
            figs[:cross_sections] = _plot_cross_sections(model)
        elseif panel === :phase_functions
            figs[:phase_functions] = _plot_phase_functions(model)
        elseif panel === :scattering
            figs[:scattering] = _plot_scattering(model)
        else
            @warn "Unknown panel: $panel. Available panels: $ALL_PANELS"
        end
    end

    return figs
end


# ======================================================================================== #
# Atmosphere
# ======================================================================================== #
function _plot_atmosphere(model)
    h_km = model.altitude_grid.h ./ 1e3
    iono = model.ionosphere
    nd = AURORA.n_neutrals(iono)

    fig = Figure(size = (800, 600))

    ax1 = Axis(fig[1, 1]; xlabel="(m⁻³)", ylabel="height (km)",
               title="MSIS neutral density", xscale=log10,
               xticks=LogTicks(LinearTicks(9)))
    lines!(ax1, nd.nN2, h_km; linewidth=2, label="N₂")
    lines!(ax1, nd.nO2, h_km; linewidth=2, linestyle=:dash, label="O₂")
    lines!(ax1, nd.nO, h_km; linewidth=2, linestyle=:dashdot, label="O")
    xlims!(ax1, 1e12, 1e20)
    ylims!(ax1, h_km[1] - 10, h_km[end] + 10)
    axislegend(ax1, "Densities"; position=:rt)

    ax2 = Axis(fig[1, 2]; xlabel="(m⁻³)", ylabel="height (km)",
               title="Electron density", xscale=log10)
    lines!(ax2, iono.ne, h_km; linewidth=2)
    xlims!(ax2, 1e8, 1e13)
    ylims!(ax2, h_km[1] - 10, h_km[end] + 10)

    return fig
end


# ======================================================================================== #
# Energy levels
# ======================================================================================== #
function _plot_energy_levels(model)
    levels = model.cross_sections.collision_levels

    fig = Figure(size = (800, 600))
    ax1 = Axis(fig[1, 1]; ylabel="Excitation energy (eV)", title="N₂")
    ax2 = Axis(fig[1, 2]; title="O₂")
    ax3 = Axis(fig[1, 3]; title="O")

    for (ax, data) in zip([ax1, ax2, ax3],
                          [levels.N2_levels, levels.O2_levels, levels.O_levels])
        for i in axes(data, 1)
            lines!(ax, [0, 1], data[i, 1] .* [1, 1]; linewidth=2)
        end
        xlims!(ax, 0, 1)
    end

    return fig
end


# ======================================================================================== #
# Energy grid
# ======================================================================================== #
function _plot_energy_grid(model)
    E = model.energy_grid.E_centers
    dE = model.energy_grid.ΔE

    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1]; xlabel="Energy (eV)", ylabel="Energy bin size (eV)",
              title="Energy variation of energy bin size")
    scatterlines!(ax, E, dE; color=:red)

    return fig
end


# ======================================================================================== #
# Cross-sections
# ======================================================================================== #
function _plot_cross_sections(model)
    E = model.energy_grid.E_centers
    σ = model.cross_sections.σ_neutrals

    fig = Figure(size=(1800, 800))

    for (idx, (species, σ_sp, title_str)) in enumerate([
            ("N2", σ.σ_N2, "σ N₂"),
            ("O2", σ.σ_O2, "σ O₂"),
            ("O",  σ.σ_O,  "σ O")])

        names = AURORA.get_level_names(species)
        n_levels = length(names)

        ax = Axis(fig[1, idx];
                    xlabel="Energy (eV)",
                    # ylabel="(m⁻²)",
                    title=title_str,
                    xscale=log10, yscale=log10,
                    xticks=LogTicks([0, 1, 2, 3, 4]),
                    xminorticksvisible=true, xminorticks=IntervalsBetween(9),
                    yticks=LogTicks(-25:-18),
                    yminorticksvisible=true, yminorticks=IntervalsBetween(9))
        if idx == 1
            ax.ylabel = "(m⁻²)"
        end
        ylims!(ax, 1e-24, 3e-19)

        # Classify each level into a visual category
        cats = [_level_category(names[i]) for i in 1:n_levels]
        n_rotvib = count(c -> c == :rotvib, cats)
        n_exc    = count(c -> c == :excitation, cats)
        n_ion    = count(c -> c == :ionization, cats)
        pal_rotvib = _make_palette(:viridis, n_rotvib)
        pal_exc    = _make_palette(:turbo, n_exc)
        pal_ion    = _make_palette(:turbo, n_ion)

        i_rv = 0; i_ex = 0; i_io = 0
        for i in 1:n_levels
            label = _format_level_name(species, names[i])
            if cats[i] == :rotvib
                i_rv += 1
                color = pal_rotvib[i_rv]
                ls = :solid; lw = 1.0
            elseif cats[i] == :ionization
                i_io += 1
                color = pal_ion[i_io]
                ls = :dash; lw = 2.0
            else
                i_ex += 1
                color = pal_exc[i_ex]
                ls = :solid; lw = 2.0
            end
            lines!(ax, E, σ_sp[i, :]; linewidth=lw, linestyle=ls, color, label)
        end

        nbanks = ceil(Int, n_levels / 4)
        Legend(fig[2, idx], ax; labelsize = 12, nbanks = nbanks,
                orientation = :horizontal,
                tellheight = true,
                tellwidth = false,
                rowgap = -2,
                valign = :top)
    end

    return fig
end

"""Format a raw level name (e.g. `"_elastic"`, `"rot0_2"`) into a readable label."""
function _format_level_name(species, name)
    # Strip leading underscore
    n = lstrip(name, '_')
    # Common substitutions for readability
    n = replace(n, "rot0_" => "rot 0→")
    n = replace(n, "vib0_" => "vib 0→")
    n = replace(n, "ionx" => "ion X")
    n = replace(n, "iona" => "ion A")
    n = replace(n, "ionb" => "ion B")
    n = replace(n, "ion16p9" => "ion 16.9eV")
    n = replace(n, "sup" => rich("Σ", subsup("u", "+")))
    n = replace(n, "sgp" => rich("Σ", subsup("g", "+")))
    n = replace(n, "sgm" => rich("Σ", subsup("g", "-")))
    n = replace(n, "pg" => rich("Π", subscript("g")))
    n = replace(n, "pu" => rich("Π", subscript("u")))
    n = replace(n, "du" => rich("Δ", subscript("u")))
    n = replace(n, "dion" => "diss. ion.")
    n = replace(n, "ddion" => "diss. diss. ion.")
    return n
end

"""
Classify a level name into `:rotvib`, `:ionization`, or `:excitation`.
"""
function _level_category(name)
    if startswith(name, "rot") || startswith(name, "vib")
        return :rotvib
    elseif startswith(name, "ion") || startswith(name, "dion") ||
           startswith(name, "ddion") || name == "dissociation"
        return :ionization
    else
        return :excitation
    end
end

"""Build a palette of `n` evenly-spaced colors from a named colormap."""
function _make_palette(cmap_name::Symbol, n)
    n <= 0 && return Makie.RGBAf[]
    cmap = Makie.cgrad(cmap_name, max(n, 2); categorical=true)
    return [cmap[(i - 1) / max(n - 1, 1)] for i in 1:n]
end


# ======================================================================================== #
# Phase functions
# ======================================================================================== #
function _plot_phase_functions(model)
    E = model.energy_grid.E_centers
    θ = deg2rad.(0:180)

    phfcnE_N2, phfcnI_N2 = AURORA.phase_fcn_N2(θ, E)
    phfcnE_O2, phfcnI_O2 = AURORA.phase_fcn_O2(θ, E)
    phfcnE_O, phfcnI_O = AURORA.phase_fcn_O(θ, E)

    custom_formatter(values) = map(
        v -> "10" * Makie.UnicodeFun.to_superscript(round(Int64, v)), values)

    fig = Figure(size=(1000, 1000))

    axis_settings = (;
        yticks=0:30:180,
        xscale=log10,
        xticks=LogTicks([1, 2, 3, 4]),
        xminorticksvisible=true,
        xminorticks=IntervalsBetween(9),
        limits=(minimum(E), maximum(E), 0, 180))
    hm_settings = (; colormap=:inferno, colorrange=(-6, -1))

    data_list = [
        (1, 1, "N₂ elastic",   phfcnE_N2),
        (1, 3, "N₂ inelastic", phfcnI_N2),
        (2, 1, "O₂ elastic",   phfcnE_O2),
        (2, 3, "O₂ inelastic", phfcnI_O2),
        (3, 1, "O elastic",    phfcnE_O),
        (3, 3, "O inelastic",  phfcnI_O),
    ]

    for (row, col, title_str, data) in data_list
        kwargs = Dict{Symbol, Any}(:title => title_str)
        if col == 1
            kwargs[:ylabel] = "Scattering angle (°)"
        end
        if row == 3
            kwargs[:xlabel] = "Energy (eV)"
        end
        ax = Axis(fig[row, col]; axis_settings..., kwargs...)
        hm = heatmap!(ax, E, 0:180, log10.(data)'; hm_settings...)
        Colorbar(fig[row, col + 1], hm; tickformat=custom_formatter,
                label="Probability")
    end

    return fig
end


# ======================================================================================== #
# Scattering matrices
# ======================================================================================== #
function _plot_scattering(model)
    P_scatter = model.scattering.P_scatter
    θ_lims = model.pitch_angle_grid.θ_lims
    n_beams = model.pitch_angle_grid.n_beams

    n_cols = ceil(Int, n_beams / 2)
    fig = Figure(size=(250 * n_cols, 500))

    for i_beam in 1:n_beams
        # Arrange: first row = first half of beams, second row = reverse of second half
        if i_beam <= n_cols
            row = 1
            col = i_beam
        else
            row = 2
            col = n_beams - i_beam + 1
        end

        θ_lo = round(min(θ_lims[i_beam], θ_lims[i_beam + 1]); digits=1)
        θ_hi = round(max(θ_lims[i_beam], θ_lims[i_beam + 1]); digits=1)

        ax = Axis(fig[row, col]; title="$(θ_hi)° – $(θ_lo)°",
                  xticks=0:30:180, yticks=0:30:180)
        if row == 2 && col == ceil(Int, n_cols / 2)
            ax.xlabel = "Scattering angle"
        end
        if col == 1
            ax.ylabel = "From angle"
        end
        heatmap!(ax, 0:180, 0:180, P_scatter[:, :, i_beam];
                 colormap=:inferno, colorrange=(0, 1))
    end

    Colorbar(fig[:, n_cols + 1]; colormap=:inferno, colorrange=(0, 1),
             label="Scattering probability")

    return fig
end
