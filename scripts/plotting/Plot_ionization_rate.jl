# This script plots the ionization rate as a function of altitude, time, and optionally
# the incoming electron flux at the top of the ionosphere.
#
# For steady-state results, `plot_excitation!` automatically switches to an altitude-profile
# line plot. For time-dependent results, a specific time step can be selected via the
# `time_index` keyword.
#
# It requires that the post-processing function `make_volume_excitation_file()` has been run
# on the simulation results first (ionization rates are saved in the same file as excitation
# rates).

using AURORA
using CairoMakie
# To interact with plots, install GLMakie and uncomment the following two lines:
# using GLMakie
# GLMakie.activate!()


## ====================================================================================== ##
## Load data
## ====================================================================================== ##

full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/backup/20260409-1812"

data = load_volume_excitation(full_path_to_directory)


## ====================================================================================== ##
## Ionization rate heatmap (time-altitude)
## ====================================================================================== ##

fig = Figure(; fontsize = 20)
ax = Axis(fig[1, 1]; xlabel = "t (s)", ylabel = "Altitude (km)",
          xminorticksvisible = true, yminorticksvisible = true,
          xticksmirrored = true, yticksmirrored = true)
hm = plot_excitation!(ax, data)
Colorbar(fig[1, 2], hm; label = "Ionization rate (/m³/s)")
display(fig)

savefile = joinpath(full_path_to_directory, "Ionization_rate_zt.png")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")


## ====================================================================================== ##
## Ionization rate + incoming flux (two-panel)
## ====================================================================================== ##

input = load_input(full_path_to_directory)

fig = Figure(; fontsize = 20)
ax_top = Axis(fig[1, 1]; yscale = log10, ylabel = "Energy (eV)",
              yminorticksvisible = true, yminorticks = IntervalsBetween(9),
              xticklabelsvisible = false, xminorticksvisible = true,
              xticksmirrored = true, yticksmirrored = true)
hm_top = plot_input!(ax_top, input, beams = 1)
Colorbar(fig[1, 2], hm_top; label = "IeE (eV/m²/s/eV/ster)")

ax_Q = Axis(fig[2, 1]; xlabel = "t (s)", ylabel = "Altitude (km)",
            xminorticksvisible = true, yminorticksvisible = true,
            xticksmirrored = true, yticksmirrored = true)
hm_Q = plot_excitation!(ax_Q, data)
Colorbar(fig[2, 2], hm_Q; label = "Ionization rate (/m³/s)")

rowsize!(fig.layout, 1, Relative(1/4))
linkxaxes!(ax_top, ax_Q)
display(fig)


## ====================================================================================== ##
## Profile mode (single time step)
## ====================================================================================== ##

## Specific time step
fig = Figure(; fontsize = 20)
ax = Axis(fig[1, 1]; xlabel = "Ionization rate (/m³/s)", ylabel = "Altitude (km)",
          xscale = log10, xminorticksvisible = true, yminorticksvisible = true,
          yminorticks = IntervalsBetween(9), xticksmirrored = true, yticksmirrored = true)
plot_excitation!(ax, data; time_index = 10)
ylims!(ax, 100, 500)
display(fig)

## Same with the input panel — mark the selected time with a vertical line
fig = Figure(; fontsize = 20)
ax_top = Axis(fig[1, 1]; yscale = log10, ylabel = "Energy (eV)",
              yminorticksvisible = true, yminorticks = IntervalsBetween(9),
              xticklabelsvisible = false, xminorticksvisible = true,
              xticksmirrored = true, yticksmirrored = true)
hm_top = plot_input!(ax_top, input)
Colorbar(fig[1, 2], hm_top; label = "IeE (eV/m²/s/eV/ster)")
vlines!(ax_top, [data.t[10]]; color = :white, linewidth = 2, linestyle = :dash)

ax_Q = Axis(fig[2, 1]; xlabel = "Ionization rate (/m³/s)", ylabel = "Altitude (km)",
            xscale = log10, xminorticksvisible = true, yminorticksvisible = true,
            yminorticks = IntervalsBetween(9), xticksmirrored = true, yticksmirrored = true)
plot_excitation!(ax_Q, data; time_index = 50)
rowsize!(fig.layout, 1, Relative(1/4))
ylims!(ax_Q, 100, 500)
display(fig)
