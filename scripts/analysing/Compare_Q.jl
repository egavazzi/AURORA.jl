using AURORA
using MAT
using CairoMakie
CairoMakie.activate!(type = "png")


directory_to_process1 = "CFL_simulations/7keV_CFL-128"
directory_to_process2 = "CFL_simulations/7keV_CFL-64"



## Load the Q data (volume emission-rates)
full_path_to_directory = pkgdir(AURORA, "data", directory_to_process1)
Q_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")
data = matread(Q_file)
Q4278_1 = data["Q4278"]
Q6730_1 = data["Q6730"]
Q7774_1 = data["Q7774"]
Q8446_1 = data["Q8446"]
QO1D_1 = data["QO1D"]
QO1S_1 = data["QO1S"]
h_atm = vec(data["h_atm"]) ./ 1e3    # convert to km
t = vec(data["t"])

## Load the Q data (volume emission-rates)
full_path_to_directory = pkgdir(AURORA, "data", directory_to_process2)
Q_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")
data = matread(Q_file)
Q4278_2 = data["Q4278"]
Q6730_2 = data["Q6730"]
Q7774_2 = data["Q7774"]
Q8446_2 = data["Q8446"]
QO1D_2 = data["QO1D"]
QO1S_2 = data["QO1S"]
h_atm = vec(data["h_atm"]) ./ 1e3    # convert to km
t = vec(data["t"])


##
diff_4278 = abs.(Q4278_1 .- Q4278_2) ./ Q4278_1
diff_6730 = abs.(Q6730_1 .- Q6730_2)
diff_7774 = abs.(Q7774_1 .- Q7774_2)
diff_8446 = abs.(Q8446_1 .- Q8446_2)

# diff_4278 = abs.(Q4278_1)
# diff_6730 = abs.(Q6730_1)
# diff_7774 = abs.(Q7774_1)
# diff_8446 = abs.(Q8446_1)

## Plot
f = with_theme(
    Theme(
        Axis = (
            xticksmirrored = true, yticksmirrored = false, xminorticksvisible = true,
            yminorticksvisible = true,
            ),
        Colorbar = (flip_vertical_label = true, vertical = true),
        Heatmap = (rasterize = true,),
        )
    ) do
    f = Figure(size = (1000, 800))
    ga = f[1, 1] = GridLayout()
    ax4278 = Axis(ga[1, 1]; title = "4278 Å", xticklabelsvisible = false, ylabel ="altitude (km)")
    hm4278 = heatmap!(t, h_atm, diff_4278'; rasterize = true)
    cb4278 = Colorbar(ga[1, 2], hm4278; label = "photons/m³/s")
    ylims!(nothing, 200)
    # xlims!(0, 0.02)
    colgap!(ga, 10)

    gb = f[1, 2] = GridLayout()
    ax6730 = Axis(gb[1, 1]; title = "6730 Å", xticklabelsvisible = false, yticklabelsvisible = false)
    hm6730 = heatmap!(t, h_atm, diff_6730')
    cb6730 = Colorbar(gb[1, 2], hm6730; label = "photons/m³/s")
    ylims!(nothing, 200)
    # xlims!(nothing, 0.05)
    colgap!(gb, 10)

    gc = f[2, 1] = GridLayout()
    ax7774 = Axis(gc[1, 1]; title = "7774 Å", xlabel = "time(s)", ylabel = "altitude (km)")
    hm7774 = heatmap!(t, h_atm, diff_7774')
    # hm7774 = contour!(t, h_atm, diff_7774'; labels=true)
    cb7774 = Colorbar(gc[1, 2], hm7774; label = "photons/m³/s")
    ylims!(nothing, 200)
    # xlims!(nothing, 0.05)
    colgap!(gc, 10)

    gd = f[2, 2] = GridLayout()
    ax8446 = Axis(gd[1, 1]; title = "8446 Å", yticklabelsvisible = false, xlabel = "time (s)")
    hm8446 = heatmap!(t, h_atm, diff_8446')
    cb8446 = Colorbar(gd[1, 2], hm8446; label = "photons/m³/s")
    ylims!(nothing, 200)
    # xlims!(nothing, 0.05)
    colgap!(gd, 10)
    return f
end

f
