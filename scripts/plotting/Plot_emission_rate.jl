# This script produces plots of:
# 1. Altitude-time variation of volume emission-rates (4278, 6730, 7774, 8446 Å)
# 2. Time variation of column-integrated emission (excitation) intensities
#
# It requires that the post-processing functions `make_volume_excitation_file()` and
# `make_column_excitation_file()` have been run on the simulation results first.

using AURORA
using CairoMakie
# To interact with plots, install GLMakie and uncomment the following two lines:
# using GLMakie
# GLMakie.activate!()


## ====================================================================================== ##
## Run
## ====================================================================================== ##

## Set the path to the simulation results directory
# full_path_to_directory = "path/to/simulation/results"
full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/backup/20260409-1812"

## Plot volume emission rates
data_Q = load_volume_excitation(full_path_to_directory)
fig_Q = plot_emission(data_Q)
display(fig_Q)

## Save volume emission plot
savefile = joinpath(full_path_to_directory, "Qtz.png")
save(savefile, fig_Q; backend = CairoMakie)
println("Saved $savefile")

## Plot column-integrated emission intensities
data_I = load_column_excitation(full_path_to_directory)
fig_I = plot_column_emission(data_I)
display(fig_I)

## Save column emission plot
savefile = joinpath(full_path_to_directory, "It.png")
save(savefile, fig_I; backend = CairoMakie)
println("Saved $savefile")


## ====================================================================================== ##
## Examples: plotting on a custom axis
## ====================================================================================== ##

## Plot a single emission wavelength on a custom axis
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "t (s)", ylabel = "Altitude (km)")
hm = plot_emission!(ax, data_Q, :Q4278)
Colorbar(fig[1, 2], hm)
fig

## Plot column emission on a custom axis
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "t (s)", ylabel = "Intensity (R)", yscale = log10)
plot_column_emission!(ax, data_I; wavelengths = [:I_4278, :I_7774])
axislegend(ax)
fig
