# This script produces plots of altitude-time-variation of volume emission-rates and
# time-variation of normalized column-emission (excitation) intensity plots.

using AURORA
using MAT
# using CairoMakie
# CairoMakie.activate!(type = "svg")
using GLMakie
GLMakie.activate!()


directory_to_plot = "Pulsating_aurora/experiment_1"



## Load the Q data (volume emission-rates)
full_path_to_directory = pkgdir(AURORA, "data", directory_to_plot)
Q_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")
data = matread(Q_file)
Q4278 = data["Q4278"]
Q6730 = data["Q6730"]
Q7774 = data["Q7774"]
Q8446 = data["Q8446"]
QO1D = data["QO1D"]
QO1S = data["QO1S"]
h_atm = vec(data["h_atm"]) ./ 1e3    # convert to km
t = vec(data["t"])

## Plot
f = with_theme(
    Theme(
        Axis = (
            xticksmirrored = true, yticksmirrored = false, xminorticksvisible = true,
            yminorticksvisible = true, limits=(nothing, (nothing, 400))
            ),
        Colorbar = (flip_vertical_label = true, vertical = true),
        Heatmap = (rasterize = true,),
        )
    ) do
    f = Figure(size = (1000, 800))
    ga = f[1, 1] = GridLayout()
    ax4278 = Axis(ga[1, 1]; title = "4278 Å", xticklabelsvisible = false, ylabel ="altitude (km)")
    hm4278 = heatmap!(t, h_atm, Q4278'; rasterize = true)
    cb4278 = Colorbar(ga[1, 2], hm4278; label = "photons/m³/s")
    colgap!(ga, 10)

    gb = f[1, 2] = GridLayout()
    ax6730 = Axis(gb[1, 1]; title = "6730 Å", xticklabelsvisible = false, yticklabelsvisible = false)
    hm6730 = heatmap!(t, h_atm, Q6730')
    cb6730 = Colorbar(gb[1, 2], hm6730; label = "photons/m³/s")
    colgap!(gb, 10)

    gc = f[2, 1] = GridLayout()
    ax7774 = Axis(gc[1, 1]; title = "7774 Å", xlabel = "time(s)", ylabel = "altitude (km)")
    hm7774 = heatmap!(t, h_atm, Q7774')
    cb7774 = Colorbar(gc[1, 2], hm7774; label = "photons/m³/s")
    colgap!(gc, 10)

    gd = f[2, 2] = GridLayout()
    ax8446 = Axis(gd[1, 1]; title = "8446 Å", yticklabelsvisible = false, xlabel = "time (s)")
    hm8446 = heatmap!(t, h_atm, Q8446')
    cb8446 = Colorbar(gd[1, 2], hm8446; label = "photons/m³/s")
    colgap!(gd, 10)
    return f
end

f

##
savefile = joinpath(full_path_to_directory, "Qtz.png")
save(savefile, f)

savefile = joinpath(full_path_to_directory, "Qtz.svg")
save(savefile, f)


##

f = Figure(size = (1000, 800))
ga = f[1, 1] = GridLayout()
ax4278 = Axis(ga[1, 1]; title = "4278 Å", xticklabelsvisible = false, ylabel ="altitude (km)")
hm4278 = heatmap!(t, h_atm, Q4278'; rasterize = true)
cb4278 = Colorbar(ga[1, 2], hm4278; label = "photons/m³/s")
ylims!(nothing, 400)
colgap!(ga, 10)
f
