# This script plots the ionization rate as a function of altitude and time, together with
# the incoming electron flux at the top of the ionosphere.
#
# It requires that the post-processing function `make_volume_excitation_file()` has been run
# on the simulation results first (ionization rates are saved in the same file as excitation
# rates).

using AURORA: beam_weight, mu_avg, find_Ietop_file
using MAT: matread
using CairoMakie
# To interact with plots, install GLMakie and uncomment the following two lines:
# using GLMakie
# GLMakie.activate!()


## ====================================================================================== ##
## Ionization rate (heatmap) + top incoming flux
## ====================================================================================== ##

function load_ionization(full_path_to_directory; plot_Ietop = true)
    # Read the ionization data
    ionization_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")
    data = matread(ionization_file)
    QO2i = data["QO2i"]
    QOi = data["QOi"]
    QN2i = data["QN2i"]
    Qtot = QN2i + QO2i + QOi
    h_atm = vec(data["h_atm"]) ./ 1e3
    t = vec(data["t"])

    # Read the Ie_top flux (optional)
    if plot_Ietop
        Ietop_file = find_Ietop_file(full_path_to_directory)
        data = matread(Ietop_file)
        Ietop = data["Ie_total"]
        μ_lims = data["mu_lims"]
        t_top = data["t_top"]; t_top = [t_top; t_top[end] + diff(t_top)[end]] .- t_top[1]
        Egrid = vec(data["E_centers"])
        dEgrid = vec(data["dE"])
    else
        Ietop = μ_lims = t_top = Egrid = dEgrid = nothing
    end

    return (; Qtot, h_atm, t, Ietop, μ_lims, t_top, Egrid, dEgrid)
end

function plot_ionization(full_path_to_directory; plot_Ietop = true)
    (; Qtot, h_atm, t, Ietop, μ_lims, t_top, Egrid, dEgrid
    ) = load_ionization(full_path_to_directory; plot_Ietop = plot_Ietop)

    fig_height = plot_Ietop ? 950 : 500
    fig = Figure(size = (850, fig_height), fontsize = 20)

    row = 1
    # Plot e- precipitation at the top of the ionosphere (optional)
    if plot_Ietop
        angle_cone = [170 180] # angle for the cone of precipitation to plot, 180° is field-aligned down
        ax_Ietop = Axis(fig[row, 1], yscale = log10, ylabel = "Energy (eV)",
                        yminorticksvisible = true, yminorticks = IntervalsBetween(9),
                        xticklabelsvisible = false, xminorticksvisible = true,
                        xticksmirrored = true, yticksmirrored = true)
        θ_lims = acosd.(μ_lims)
        idx_θ = vec(angle_cone[1] .<= abs.(acosd.(mu_avg(θ_lims))) .<= angle_cone[2])
        BeamW = beam_weight([angle_cone[1], angle_cone[2]])
        data_heatmap = dropdims(sum(Ietop[idx_θ, :, :]; dims=1); dims=1) ./ BeamW ./ dEgrid' .* Egrid'
        hm_Ietop = heatmap!(t_top, Egrid, data_heatmap;
                            colorrange = (1e6, maximum(data_heatmap)),
                            colorscale = log10, colormap = :inferno)
        Colorbar(fig[row, 2], hm_Ietop; label = "IeE (eV/m²/s/eV/ster)")
        row += 1
    end

    # Plot ionization rate as heatmap
    ax_Q = Axis(fig[row, 1], xlabel = "t (s)", ylabel = "Altitude (km)",
                xminorticksvisible = true, yminorticksvisible = true,
                xticksmirrored = true, yticksmirrored = true)
    hm_Q = heatmap!(t, h_atm, Qtot'; colorrange = (1e6, maximum(Qtot)),
                    colorscale = log10, rasterize = true)
    Colorbar(fig[row, 2], hm_Q; label = "Ionization rate (/m³/s)")
    # Contour overlay
    Qtot_contour = log10.(Qtot) |> (x -> replace(x, -Inf => 0))
    contour!(t, h_atm, Qtot_contour'; levels = 6:10, color = :black, labels = true)

    # Layout adjustments
    colsize!(fig.layout, 1, Auto())
    if plot_Ietop
        rowsize!(fig.layout, 1, Relative(1/5))
        linkxaxes!(ax_Ietop, ax_Q)
    end
    xlims!(ax_Q, 0, min(plot_Ietop ? t_top[end] : t[end], t[end]))
    ylims!(ax_Q, 100, 500)

    # Height of maximum ionization
    height_Qmax = [h_atm[i_max[1]] for i_max in vec(findmax(Qtot, dims=1)[2])]
    height_Qmax[vec(maximum(Qtot, dims=1)) .< maximum(Qtot) / 10] .= NaN
    lines!(ax_Q, t, height_Qmax; color = :red, linestyle = :dash, linewidth = 2)

    return fig
end


## ====================================================================================== ##
## Run
## ====================================================================================== ##

## Set the path to the simulation results directory
# full_path_to_directory = "path/to/simulation/results"
full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/backup/20260409-1812"

## Plot ionization rate
fig = plot_ionization(full_path_to_directory, plot_Ietop = false)
display(fig)

## Save the figure
savefile = joinpath(full_path_to_directory, "Ionization_rate_zt.png")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")
