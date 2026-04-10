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
## Load data
## ====================================================================================== ##

full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/backup/20260410-1015"
full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/backup/20260409-1812"

data_Q = load_volume_excitation(full_path_to_directory)
data_I = load_column_excitation(full_path_to_directory)


## ====================================================================================== ##
## Volume emission rates (2×2 heatmap grid)
## ====================================================================================== ##

wavelengths = [:Q4278, :Q6730, :Q7774, :Q8446]
positions   = [(1, 1), (1, 2), (2, 1), (2, 2)]
ax_labels   = Dict(:Q4278 => "4278 Å", :Q6730 => "6730 Å",
                   :Q7774 => "7774 Å", :Q8446 => "8446 Å")

fig = Figure(; fontsize = 20)
ax_kwargs = (; xticksmirrored = true, yticksmirrored = true,
               xminorticksvisible = true, yminorticksvisible = true,
               limits = ((0, data_Q.t[end]), (100, 400)))

for (wl, (r, c)) in zip(wavelengths, positions)
    g = fig[r, c] = GridLayout()
    ax = Axis(g[1, 1]; ax_kwargs..., title = ax_labels[wl],
              xlabel = r == 2 ? "t (s)" : "",
              ylabel = c == 1 ? "Altitude (km)" : "",
              xticklabelsvisible = r == 2,
              yticklabelsvisible = c == 1)
    hm = plot_excitation!(ax, data_Q; field = wl)
    Colorbar(g[1, 2], hm; label = "photons/m³/s")
    colgap!(g, 10)
end
display(fig)

savefile = joinpath(full_path_to_directory, "Qtz.png")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")


## ====================================================================================== ##
## Column-integrated emission intensities
## ====================================================================================== ##

fig = Figure(; fontsize = 20)
ax = Axis(fig[1, 1]; xlabel = "t (s)", ylabel = "Intensity (R)",
          yscale = log10,
          xminorticksvisible = true, xminorgridvisible = true,
          xminorticks = IntervalsBetween(10),
          yminorticksvisible = true, yminorgridvisible = true,
          yminorticks = IntervalsBetween(9), yticksmirrored = true)
plot_column_excitation!(ax, data_I)
Legend(fig[1, 2], ax; patchsize = [40, 20])
display(fig)

savefile = joinpath(full_path_to_directory, "It.png")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")


## ====================================================================================== ##
## Profile mode (single time step)
## ====================================================================================== ##

fig = Figure(; fontsize = 20)
ax_kwargs_profile = (; xscale = log10, xticksmirrored = true, yticksmirrored = true,
                       xminorticksvisible = true, yminorticksvisible = true,
                       yminorticks = IntervalsBetween(9),
                       limits = (nothing, (100, 400)))
for (wl, (r, c)) in zip(wavelengths, positions)
    ax = Axis(fig[r, c]; ax_kwargs_profile..., title = ax_labels[wl],
              xlabel = r == 2 ? "photons/m³/s" : "",
              ylabel = c == 1 ? "Altitude (km)" : "",
              xticklabelsvisible = r == 2,
              yticklabelsvisible = c == 1)
    plot_excitation!(ax, data_Q; field = wl, time_index = 1)
end
display(fig)
