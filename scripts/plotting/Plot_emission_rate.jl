# This script produces plots of:
# 1. Altitude-time variation of volume emission-rates (4278, 6730, 7774, 8446 Å)
# 2. Time variation of column-integrated emission (excitation) intensities
#
# It requires that the post-processing functions `make_volume_excitation_file()` and
# `make_column_excitation_file()` have been run on the simulation results first.

using AURORA: beam_weight, mu_avg
using MAT: matread
using CairoMakie
# To interact with plots, install GLMakie and uncomment the following two lines:
# using GLMakie
# GLMakie.activate!()


## ====================================================================================== ##
## Volume emission rates (altitude × time heatmaps)
## ====================================================================================== ##

function load_volume_emission(full_path_to_directory)
    Q_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")
    data = matread(Q_file)
    Q4278 = data["Q4278"]
    Q6730 = data["Q6730"]
    Q7774 = data["Q7774"]
    Q8446 = data["Q8446"]
    h_atm = vec(data["h_atm"]) ./ 1e3 # convert to km
    t = vec(data["t"])

    # Find height of maximum emission for each wavelength
    function height_of_max(Q, h)
        h_max = [h[i_max[1]] for i_max in vec(findmax(Q, dims=1)[2])]
        h_max[vec(maximum(Q, dims=1)) .< maximum(Q) / 10] .= NaN
        return h_max
    end

    return (; Q4278, Q6730, Q7774, Q8446, h_atm, t,
            height_Q4278max = height_of_max(Q4278, h_atm),
            height_Q6730max = height_of_max(Q6730, h_atm),
            height_Q7774max = height_of_max(Q7774, h_atm),
            height_Q8446max = height_of_max(Q8446, h_atm))
end

function plot_volume_emission(full_path_to_directory)
    (; Q4278, Q6730, Q7774, Q8446, h_atm, t,
     height_Q4278max, height_Q6730max, height_Q7774max, height_Q8446max
    ) = load_volume_emission(full_path_to_directory)

    fig = Figure(size = (1000, 800), fontsize = 20)

    ax_kwargs = (; xticksmirrored = true, yticksmirrored = true,
                   xminorticksvisible = true, yminorticksvisible = true,
                   limits = ((0, t[end]), (100, 400)))

    colorrange = (1e6, 1e7) # log10 of the color range for the heatmaps
    ga = fig[1, 1] = GridLayout()
    ax4278 = Axis(ga[1, 1]; ax_kwargs..., title = "4278 Å",
                  xticklabelsvisible = false, ylabel = "Altitude (km)")
    hm4278 = heatmap!(t, h_atm, Q4278'; colorscale = log10, colorrange)
    Colorbar(ga[1, 2], hm4278; label = "photons/m³/s")
    lines!(ax4278, t, height_Q4278max; color = :red, linestyle = :dash, linewidth = 2)
    colgap!(ga, 10)

    gb = fig[1, 2] = GridLayout()
    ax6730 = Axis(gb[1, 1]; ax_kwargs..., title = "6730 Å",
                  xticklabelsvisible = false, yticklabelsvisible = false)
    hm6730 = heatmap!(t, h_atm, Q6730'; colorscale = log10, colorrange)
    Colorbar(gb[1, 2], hm6730; label = "photons/m³/s")
    lines!(ax6730, t, height_Q6730max; color = :red, linestyle = :dash, linewidth = 2)
    colgap!(gb, 10)

    gc = fig[2, 1] = GridLayout()
    ax7774 = Axis(gc[1, 1]; ax_kwargs..., title = "7774 Å",
                  xlabel = "t (s)", ylabel = "Altitude (km)")
    hm7774 = heatmap!(t, h_atm, Q7774'; colorscale = log10, colorrange)
    Colorbar(gc[1, 2], hm7774; label = "photons/m³/s")
    lines!(ax7774, t, height_Q7774max; color = :red, linestyle = :dash, linewidth = 2)
    colgap!(gc, 10)

    gd = fig[2, 2] = GridLayout()
    ax8446 = Axis(gd[1, 1]; ax_kwargs..., title = "8446 Å",
                  yticklabelsvisible = false, xlabel = "t (s)")
    hm8446 = heatmap!(t, h_atm, Q8446'; colorscale = log10, colorrange)
    Colorbar(gd[1, 2], hm8446; label = "photons/m³/s")
    lines!(ax8446, t, height_Q8446max; color = :red, linestyle = :dash, linewidth = 2)
    colgap!(gd, 10)

    return fig
end


## ====================================================================================== ##
## Column-integrated emission intensities (line plot)
## ====================================================================================== ##

function load_column_emission(full_path_to_directory)
    I_file = joinpath(full_path_to_directory, "I_lambda_of_t.mat")
    data = matread(I_file)
    # Divide by 1e10 to convert from photons/m² to Rayleigh
    I4278 = vec(data["I_4278"]) / 1e10
    I6730 = vec(data["I_6730"]) / 1e10
    I7774 = vec(data["I_7774"]) / 1e10
    I8446 = vec(data["I_8446"]) / 1e10
    IO1D  = vec(data["I_O1D"]) / 1e10
    IO1S  = vec(data["I_O1S"]) / 1e10
    t = vec(data["t"])
    return (; I4278, I6730, I7774, I8446, IO1D, IO1S, t)
end

function plot_column_emission(full_path_to_directory)
    (; I4278, I6730, I7774, I8446, IO1D, IO1S, t
    ) = load_column_emission(full_path_to_directory)

    fig = Figure(size = (1000, 800), fontsize = 20)
    ax = Axis(fig[1, 1]; xlabel = "t (s)", ylabel = "Intensity (R)",
              yscale = log10,
              xminorticksvisible = true, xminorgridvisible = true,
              xminorticks = IntervalsBetween(10),
              yminorticksvisible = true, yminorgridvisible = true,
              yminorticks = IntervalsBetween(9), yticksmirrored = true)
    lines!(t, I4278; label = rich("I", subscript("4278")), color = :blue)
    lines!(t, I6730; label = rich("I", subscript("6730")), color = :red)
    lines!(t, I7774; label = rich("I", subscript("7774")), color = RGBf(0.5, 0, 0))
    lines!(t, I8446; label = rich("I", subscript("8446")), color = :black)
    lines!(t, IO1D;  label = rich("q", subscript("O(¹D)")), color = RGBf(1, 0.2, 0), linestyle = :dash)
    lines!(t, IO1S;  label = rich("q", subscript("O(¹S)")), color = :green, linestyle = :dash)
    Legend(fig[1, 2], ax; patchsize = [40, 20])

    return fig
end


## ====================================================================================== ##
## Run
## ====================================================================================== ##

## Set the path to the simulation results directory
# full_path_to_directory = "path/to/simulation/results"
full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/backup/20260409-1812"

## Plot volume emission rates
fig_Q = plot_volume_emission(full_path_to_directory)
display(fig_Q)

## Save volume emission plot
savefile = joinpath(full_path_to_directory, "Qtz.png")
save(savefile, fig_Q; backend = CairoMakie)
println("Saved $savefile")

## Plot column-integrated emission intensities
fig_I = plot_column_emission(full_path_to_directory)
display(fig_I)

## Save column emission plot
savefile = joinpath(full_path_to_directory, "It.png")
save(savefile, fig_I; backend = CairoMakie)
println("Saved $savefile")
