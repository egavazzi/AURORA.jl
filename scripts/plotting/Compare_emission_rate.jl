using AURORA
using MAT
using CairoMakie
CairoMakie.activate!()
using GLMakie
# GLMakie.activate!()


## This is for simple cases of comparison where z-grid and t-grid are the same
# Path to dta folder 1
full_path_to_directory1 = joinpath(REVONTULI_MOUNT, "mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_536s_lr_Bz-9_newZ_550km_finer-theta_halfstepsAB")
# Path to data folder 2
full_path_to_directory2 = joinpath(REVONTULI_MOUNT, "mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_536s_lr_Bz-9_newZ_550km_finer-theta_halfstepsAB_scaled")


## Read the Q data (volume emission-rates)
# From data folder 1
Q_file = joinpath(full_path_to_directory1, "Qzt_all_L.mat")
data = matread(Q_file)
Q4278_1 = data["Q4278"]
Q6730_1 = data["Q6730"]
Q7774_1 = data["Q7774"]
Q8446_1 = data["Q8446"]
QO1D_1 = data["QO1D"]
QO1S_1 = data["QO1S"]
h_atm = vec(data["h_atm"]) ./ 1e3 # convert to km
t = vec(data["t"])
# From data folder 2
Q_file = joinpath(full_path_to_directory2, "Qzt_all_L.mat")
data = matread(Q_file)
Q4278_2 = data["Q4278"]
Q6730_2 = data["Q6730"]
Q7774_2 = data["Q7774"]
Q8446_2 = data["Q8446"]
QO1D_2 = data["QO1D"]
QO1S_2 = data["QO1S"]
# h_atm and t are not loaded again (they are the same)


## Plot Q data + difference
fig = with_theme(
    Theme(
        Axis = (
            xticksmirrored = true, yticksmirrored = true, xminorticksvisible = true,
            yminorticksvisible = true, limits=(nothing, (nothing, 400))
            ),
            Heatmap = (rasterize = true,),
            Colorbar = (flip_vertical_label = true, vertical = true),
        )
    ) do

    fig = Figure(size = (1400, 1200))
    Label(fig[0, 1:2],
          "1st column: " * split(full_path_to_directory1, "/")[end] * "\n" *
          "2nd column: " * split(full_path_to_directory2, "/")[end],
          tellheight = true, tellwidth = false, font = :bold, halign = :left,
          justification = :left)

    # First row : 4278 Å
    colorlims = (0, maximum([maximum(Q4278_1), maximum(Q4278_2)]))
    colorlims_diff = (-maximum(abs.(Q4278_1 - Q4278_2)), maximum(abs.(Q4278_1 - Q4278_2)))
    ga_1 = fig[1, 1] = GridLayout()
    ax4278 = Axis(ga_1[1, 1]; title = "4278 Å",
                  xticklabelsvisible = false, ylabel = "altitude (km)")
    hm4278 = heatmap!(t, h_atm, Q4278_1'; colorrange = colorlims)
    cb4278 = Colorbar(ga_1[1, 2], hm4278; label = "photons/m³/s")
    colgap!(ga_1, 10)
    ga_2 = fig[1, 2] = GridLayout()
    ax4278 = Axis(ga_2[1, 1]; title = "4278 Å",
                  xticklabelsvisible = false, yticklabelsvisible = false)
    hm4278 = heatmap!(t, h_atm, Q4278_2'; colorrange = colorlims)
    cb4278 = Colorbar(ga_2[1, 2], hm4278; label = "photons/m³/s")
    colgap!(ga_2, 10)
    ga_3 = fig[1, 3] = GridLayout()
    ax4278 = Axis(ga_3[1, 1]; title = "Difference (1 - 2)",
                  xticklabelsvisible = false, yticklabelsvisible = false)
    hm4278 = heatmap!(t, h_atm, (Q4278_1 - Q4278_2)';
                      colormap = :RdBu, colorrange = colorlims_diff)
    cb4278 = Colorbar(ga_3[1, 2], hm4278; label = "photons/m³/s")
    colgap!(ga_3, 10)
    ga_4 = fig[1, 4] = GridLayout()
    ax4278 = Axis(ga_4[1, 1]; title = "Relative difference (1 - 2)",
    xticklabelsvisible = false, yticklabelsvisible = false)
    hm4278 = heatmap!(t, h_atm, (abs.(Q4278_1 - Q4278_2) ./ Q4278_1)';
            colormap = :cividis, colorrange = (0, 1))
    cb4278 = Colorbar(ga_4[1, 2], hm4278; label = "relative difference")
    colgap!(ga_4, 10)

    # Second row : 6730 Å
    colorlims = (0, maximum([maximum(Q6730_1), maximum(Q6730_2)]))
    colorlims_diff = (-maximum(abs.(Q6730_1 - Q6730_2)), maximum(abs.(Q6730_1 - Q6730_2)))
    gb_1 = fig[2, 1] = GridLayout()
    ax6730 = Axis(gb_1[1, 1]; title = "6730 Å",
                  xticklabelsvisible = false, ylabel = "altitude (km)")
    hm6730 = heatmap!(t, h_atm, Q6730_1'; colorrange = colorlims)
    cb6730 = Colorbar(gb_1[1, 2], hm6730; label = "photons/m³/s")
    colgap!(gb_1, 10)
    gb_2 = fig[2, 2] = GridLayout()
    ax6730 = Axis(gb_2[1, 1]; title = "6730 Å",
                  xticklabelsvisible = false, yticklabelsvisible = false)
    hm6730 = heatmap!(t, h_atm, Q6730_2'; colorrange = colorlims)
    cb6730 = Colorbar(gb_2[1, 2], hm6730; label = "photons/m³/s")
    colgap!(gb_2, 10)
    gb_3 = fig[2, 3] = GridLayout()
    ax6730 = Axis(gb_3[1, 1]; title = "Difference (1 - 2)",
                  xticklabelsvisible = false, yticklabelsvisible = false)
    hm6730 = heatmap!(t, h_atm, (Q6730_1 - Q6730_2)';
                      colormap = :RdBu, colorrange = colorlims_diff)
    cb6730 = Colorbar(gb_3[1, 2], hm6730; label = "photons/m³/s")
    colgap!(gb_3, 10)
    gb_4 = fig[2, 4] = GridLayout()
    ax6730 = Axis(gb_4[1, 1]; title = "Relative difference (1 - 2)",
    xticklabelsvisible = false, yticklabelsvisible = false)
    hm6730 = heatmap!(t, h_atm, (abs.(Q6730_1 - Q6730_2) ./ Q6730_1)';
            colormap = :cividis, colorrange = (0, 1))
    cb6730 = Colorbar(gb_4[1, 2], hm6730; label = "relative difference")
    colgap!(gb_4, 10)


    # Third row : 7774 Å
    colorlims = (0, maximum([maximum(Q7774_1), maximum(Q7774_2)]))
    colorlims_diff = (-maximum(abs.(Q7774_1 - Q7774_2)), maximum(abs.(Q7774_1 - Q7774_2)))
    gc_1 = fig[3, 1] = GridLayout()
    ax7774 = Axis(gc_1[1, 1]; title = "7774 Å",
                  xticklabelsvisible = false, ylabel = "altitude (km)")
    hm7774 = heatmap!(t, h_atm, Q7774_1'; colorrange = colorlims)
    cb7774 = Colorbar(gc_1[1, 2], hm7774; label = "photons/m³/s")
    colgap!(gc_1, 10)
    gc_2 = fig[3, 2] = GridLayout()
    ax7774 = Axis(gc_2[1, 1]; title = "7774 Å",
                  xticklabelsvisible = false, yticklabelsvisible = false)
    hm7774 = heatmap!(t, h_atm, Q7774_2'; colorrange = colorlims)
    cb7774 = Colorbar(gc_2[1, 2], hm7774; label = "photons/m³/s")
    colgap!(gc_2, 10)
    gc_3 = fig[3, 3] = GridLayout()
    ax7774 = Axis(gc_3[1, 1]; title = "Difference (1 - 2)",
                  xticklabelsvisible = false, yticklabelsvisible = false)
    hm7774 = heatmap!(t, h_atm, (Q7774_1 - Q7774_2)';
                      colormap = :RdBu, colorrange = colorlims_diff)
    cb7774 = Colorbar(gc_3[1, 2], hm7774; label = "photons/m³/s")
    colgap!(gc_3, 10)
    gc_4 = fig[3, 4] = GridLayout()
    ax7774 = Axis(gc_4[1, 1]; title = "Relative difference (1 - 2)",
    xticklabelsvisible = false, yticklabelsvisible = false)
    hm7774 = heatmap!(t, h_atm, (abs.(Q7774_1 - Q7774_2) ./ Q7774_1)';
            colormap = :cividis, colorrange = (0, 1))
    cb7774 = Colorbar(gc_4[1, 2], hm7774; label = "relative difference")
    colgap!(gc_4, 10)

    # Fourth row : 8446 Å
    colorlims = (0, maximum([maximum(Q8446_1), maximum(Q8446_2)]))
    colorlims_diff = (-maximum(abs.(Q8446_1 - Q8446_2)), maximum(abs.(Q8446_1 - Q8446_2)))
    gd_1 = fig[4, 1] = GridLayout()
    ax8446 = Axis(gd_1[1, 1]; title = "8446 Å",
                  xlabel = "time (s)", ylabel = "altitude (km)")
    hm8446 = heatmap!(t, h_atm, Q8446_1'; colorrange = colorlims)
    cb8446 = Colorbar(gd_1[1, 2], hm8446; label = "photons/m³/s")
    colgap!(gd_1, 10)
    gd_2 = fig[4, 2] = GridLayout()
    ax8446 = Axis(gd_2[1, 1]; title = "8446 Å",
                  xlabel = "time (s)", yticklabelsvisible = false)
    hm8446 = heatmap!(t, h_atm, Q8446_2'; colorrange = colorlims)
    cb8446 = Colorbar(gd_2[1, 2], hm8446; label = "photons/m³/s")
    colgap!(gd_2, 10)
    gd_3 = fig[4, 3] = GridLayout()
    ax8446 = Axis(gd_3[1, 1]; title = "Difference (1 - 2)",
                  xlabel = "time (s)", yticklabelsvisible = false)
    hm8446 = heatmap!(t, h_atm, (Q8446_1 - Q8446_2)';
                      colormap = :RdBu, colorrange = colorlims_diff)
    cb8446 = Colorbar(gd_3[1, 2], hm8446; label = "photons/m³/s")
    colgap!(gd_3, 10)
    gd_4 = fig[4, 4] = GridLayout()
    ax8446 = Axis(gd_4[1, 1]; title = "Relative difference (1 - 2)",
    xticklabelsvisible = false, yticklabelsvisible = false)
    hm8446 = heatmap!(t, h_atm, (abs.(Q8446_1 - Q8446_2) ./ Q8446_1)';
            colormap = :cividis, colorrange = (0, 1))
    cb8446 = Colorbar(gd_4[1, 2], hm8446; label = "relative difference")
    colgap!(gd_4, 10)

    return fig
end
display(fig)
# display(GLMakie.Screen(), fig)
