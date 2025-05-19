# Kind of a WIP

using AURORA
using MAT
using CairoMakie
CairoMakie.activate!()
# using GLMakie
# GLMakie.activate!()


## This is for simple cases of comparison where z-grid and t-grid are the same
# Path to data folder 1
# full_path_to_data1 = "data/backup/20250504-2044/IeFlickering-01.mat"
full_path_to_data1 = "data/compare/new_1keV_constant_CFL256/IeFlickering-01.mat"
# Path to data folder 2
# full_path_to_data2 = "data/backup/20250504-2027/IeFlickering-01.mat"
full_path_to_data2 = "data/compare/new_1keV_oldQ_constant_CFL256/IeFlickering-01.mat"

# Read the data (electron number flux)
# From data folder 1
data = matread(full_path_to_data1)
Ie1 = data["Ie_ztE"]
h_atm = vec(data["h_atm"]) ./ 1e3 # convert to km
t = vec(data["t_run"])
E = data["E"]
μ_lims = data["mu_lims"]
# From data folder 2
data = matread(full_path_to_data2)
Ie2 = data["Ie_ztE"]


Ie1_4D = AURORA.restructure_Ie_from_3D_to_4D(Ie1, μ_lims, h_atm, t, E);
Ie2_4D = AURORA.restructure_Ie_from_3D_to_4D(Ie2, μ_lims, h_atm, t, E);

# Plot Ie data + difference (fixed time, in height and energy)
i_μ = 1
i_t = length(t)

data_plot1 = Ie1_4D[i_μ, :, i_t, :];
data_plot2 = Ie2_4D[i_μ, :, i_t, :];

##
Theme(Axis = (xticksmirrored = true, yticksmirrored = true, xminorticksvisible = true,
              yminorticksvisible = true, limits = (nothing, (nothing, 400))),
      Heatmap = (rasterize = true,),
      Colorbar = (flip_vertical_label = true, vertical = true))

fig = Figure(size = (1600, 800))
Label(fig[0, 1:2],
        "Ie1: " * split(full_path_to_data1, "/")[end - 1] * "\n" *
        "Ie2: " * split(full_path_to_data2, "/")[end - 1],
        tellheight = true, tellwidth = false, font = :bold, halign = :left,
        justification = :left)

# First row : 4278 Å
colorlims = (0, maximum([maximum(data_plot1), maximum(data_plot2)]))
# colorlims_diff = (-maximum(abs.(data_plot1 - data_plot2)), maximum(abs.(data_plot1 - data_plot2)))
colorlims_diff = (1e-10, maximum(abs.(data_plot1 - data_plot2)))
ga_1 = fig[1, 1] = GridLayout()
ax1 = Axis(ga_1[1, 1]; title = "Ie1", ylabel = "altitude (km)")
hm1 = heatmap!(E, h_atm, data_plot1'; colorrange = colorlims)
cb1 = Colorbar(ga_1[1, 2], hm1; label = "")
colgap!(ga_1, 10)
ga_2 = fig[2, 1] = GridLayout()
ax2 = Axis(ga_2[1, 1]; title = "Ie2", ylabel = "altitude (km)")
hm2 = heatmap!(E, h_atm, data_plot2'; colorrange = colorlims)
cb2 = Colorbar(ga_2[1, 2], hm2; label = "")
colgap!(ga_2, 10)
if data_plot1 != data_plot2
    ga_3 = fig[1:2, 2:3] = GridLayout()
    ax3 = Axis(ga_3[1, 1]; title = "Difference (1 - 2)",
                    yticklabelsvisible = false)
    hm3 = heatmap!(E, h_atm, abs.(data_plot1 - data_plot2)';
                        colormap = :RdBu, colorscale=log10, colorrange = colorlims_diff)
    cb3 = Colorbar(ga_3[1, 2], hm3; label = "")
    colgap!(ga_3, 10)
    ga_4 = fig[1:2, 4] = GridLayout()
    ax4 = Axis(ga_4[1, 1]; title = "Relative difference (1 - 2)",
                    yticklabelsvisible = false)
    hm4 = heatmap!(E, h_atm, (abs.(data_plot1 - data_plot2) ./ data_plot1)';
            colormap = :cividis, colorrange = (0, 1))
    cb4 = Colorbar(ga_4[1, 2], hm4; label = "relative difference")
    colgap!(ga_4, 10)

    # xlims!(ax1, 0, 20)
    # xlims!(ax2, 0, 20)
    # xlims!(ax3, 0, 20)
else
    Label(fig[:, 2], "The results are exactly the same", tellwidth = false)
end

display(fig)
# display(GLMakie.Screen(), fig)





##

# Plot Ie data + difference (fixed energy, in height and time)
i_μ = 1
iE = 10
data_plot1 = Ie1_4D[i_μ, :, :, iE];
data_plot2 = Ie2_4D[i_μ, :, :, iE];

#
Theme(Axis = (xticksmirrored = true, yticksmirrored = true, xminorticksvisible = true,
              yminorticksvisible = true, limits = (nothing, (nothing, 400))),
      Heatmap = (rasterize = true,),
      Colorbar = (flip_vertical_label = true, vertical = true))

fig = Figure(size = (1600, 800))
Label(fig[0, 1:2],
        "Ie1: " * split(full_path_to_data1, "/")[end - 1] * "\n" *
        "Ie2: " * split(full_path_to_data2, "/")[end - 1],
        tellheight = true, tellwidth = false, font = :bold, halign = :left,
        justification = :left)

Label(fig[0, 4],
        "E = $(E[iE])",
        tellheight = true, tellwidth = false, font = :bold, halign = :left,
        justification = :left)

# First row : 4278 Å
colorlims = (0, maximum([maximum(data_plot1), maximum(data_plot2)]))
# colorlims_diff = (-maximum(abs.(data_plot1 - data_plot2)), maximum(abs.(data_plot1 - data_plot2)))
colorlims_diff = (1e-6, maximum(abs.(data_plot1 - data_plot2)))
ga_1 = fig[1, 1] = GridLayout()
ax1 = Axis(ga_1[1, 1]; title = "Ie1", ylabel = "altitude (km)")
hm1 = heatmap!(t, h_atm, data_plot1'; colorrange = colorlims)
cb1 = Colorbar(ga_1[1, 2], hm1; label = "")
colgap!(ga_1, 10)
ga_2 = fig[2, 1] = GridLayout()
ax2 = Axis(ga_2[1, 1]; title = "Ie2", ylabel = "altitude (km)")
hm2 = heatmap!(t, h_atm, data_plot2'; colorrange = colorlims)
cb2 = Colorbar(ga_2[1, 2], hm2; label = "")
colgap!(ga_2, 10)
if data_plot1 != data_plot2
    ga_3 = fig[1:2, 2:3] = GridLayout()
    ax3 = Axis(ga_3[1, 1]; title = "Difference (1 - 2)",
                    yticklabelsvisible = false)
    hm3 = heatmap!(t, h_atm, abs.(data_plot1 - data_plot2)';
                        colormap = :RdBu, colorscale=log10, colorrange = colorlims_diff)
    cb3 = Colorbar(ga_3[1, 2], hm3; label = "")
    colgap!(ga_3, 10)
    ga_4 = fig[1:2, 4] = GridLayout()
    ax4 = Axis(ga_4[1, 1]; title = "Relative difference (1 - 2)",
                    yticklabelsvisible = false)
    hm4 = heatmap!(t, h_atm, (abs.(data_plot1 - data_plot2) ./ data_plot1)';
            colormap = :cividis, colorrange = (0, 1))
    cb4 = Colorbar(ga_4[1, 2], hm4; label = "relative difference")
    colgap!(ga_4, 10)

    # xlims!(ax1, 0, 20)
    # xlims!(ax2, 0, 20)
    # xlims!(ax3, 0, 20)
else
    Label(fig[:, 2], "The results are exactly the same", tellwidth = false)
end

display(fig)
# display(GLMakie.Screen(), fig)
