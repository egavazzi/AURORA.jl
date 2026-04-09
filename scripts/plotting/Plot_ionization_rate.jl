# This script plots the ionization rate as a function of altitude and time, together with
# the incoming electron flux at the top of the ionosphere.
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
## Run
## ====================================================================================== ##

## Set the path to the simulation results directory
# full_path_to_directory = "path/to/simulation/results"
full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/backup/20260409-1812"

## Plot ionization rate
data = load_volume_excitation(full_path_to_directory)
fig = plot_ionization(data; plot_input = true)
display(fig)

## Save the figure
savefile = joinpath(full_path_to_directory, "Ionization_rate_zt.png")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")


## ====================================================================================== ##
## Examples: plotting on a custom axis
## ====================================================================================== ##

## Plot ionization rate on a custom axis (e.g. for a single species)
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "t (s)", ylabel = "Altitude (km)")
hm = plot_ionization!(ax, data; species = :total)
Colorbar(fig[1, 2], hm)
fig
