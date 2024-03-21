# This script produces plots of altitude-time-variation of volume emission-rates and
# time-variation of normalized column-emission (excitation) intensity plots.

using AURORA
using MAT
using GLMakie
# GLMakie.activate!()
using CairoMakie
CairoMakie.activate!()


# directory to plot, absolute path
# full_path_to_directory = joinpath(REVONTULI_MOUNT, "mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_536s_lr_Bz-9_newZ_550km_finer-theta_halfstepsAB_scaled")
full_path_to_directory = joinpath(REVONTULI_MOUNT, "mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_536s_correct_msis_and_scattering")



## Load the Q data (volume emission-rates)
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



## Plot Q data
f = with_theme(
    Theme(
        Axis = (
            xticksmirrored = true, yticksmirrored = false, xminorticksvisible = true,
            yminorticksvisible = true, limits=(nothing, (nothing, 400))
            ),
            Heatmap = (rasterize = true,),
            Colorbar = (flip_vertical_label = true, vertical = true),
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
display(f)
# display(GLMakie.Screen(), f)



## save plot as image
savefile = joinpath(full_path_to_directory, "Qtz.png")
save(savefile, f)

savefile = joinpath(full_path_to_directory, "Qtz.svg")
save(savefile, f)
##














## Load the I data (column-integrated emission-rates)
I_file = joinpath(full_path_to_directory, "I_lambda_of_t.mat")
data = matread(I_file)
I4278 = vec(data["I_4278"])
I6730 = vec(data["I_6730"])
I7774 = vec(data["I_7774"])
I7774_O = vec(data["I_7774_O"])
I7774_O2 = vec(data["I_7774_O2"])
I8446 = vec(data["I_8446"])
I8446_O = vec(data["I_8446_O"])
I8446_O2 = vec(data["I_8446_O2"])
IO1D = vec(data["I_O1D"])
IO1S = vec(data["I_O1S"])
t = vec(data["t"])



## Plot the I data
custom_theme = Theme(fontsize = 20, linewidth = 2)
set_theme!(custom_theme)
f = Figure(size = (1000, 800))
ax = Axis(f[1, 1]; title = "Intensity", xlabel = "time (s)", ylabel = "#exc/m²/s",
          yscale = log10, limits = (t[1], t[end], 1e4, 1e11), xticks = t[1]:0.1:t[end],
          yminorticksvisible = true, yminorgridvisible = true,
          yminorticks = IntervalsBetween(9), yticksmirrored = true)
ylims!(1e4, 1e11)
lines!(t[2:end], I4278[2:end]; label = rich("I", subscript("4278")), color = :blue)
lines!(t[2:end], I6730[2:end]; label = rich("I", subscript("6730")), color = :red)
lines!(t, I7774; label = rich("I", subscript("7774")), color = RGBf(0.5, 0, 0))
# lines!(t, I7774_O; label = rich("I", subscript("7774(O)")))
# lines!(t, I7774_O2; label = rich("I", subscript("7774(O2)")))
lines!(t, I8446; label = rich("I", subscript("8446")), color = :black)
# lines!(t, I8446_O; label = rich("I", subscript("8446(O)")))
# lines!(t, I8446_O2; label = rich("I", subscript("8446(O2)")))
lines!(t, IO1D; label = rich("I", subscript("O(¹D)")), color = RGBf(1, 0.2, 0), linestyle = :dash)
lines!(t, IO1S; label = rich("I", subscript("O(¹S)")), color = :green, linestyle = :dash)
# axislegend(ax, position = :lt)
Legend(f[1, 2], ax; patchsize = [40, 20])
set_theme!()
display(f)
# display(GLMakie.Screen(), f)



## save plot as image
savefile = joinpath(full_path_to_directory, "It.png")
save(savefile, f)

savefile = joinpath(full_path_to_directory, "It.svg")
save(savefile, f)
##




























##

f = Figure(size = (1000, 800))
ga = f[1, 1] = GridLayout()
ax4278 = Axis(ga[1, 1]; title = "4278 Å", xticklabelsvisible = false, ylabel ="altitude (km)")
hm4278 = heatmap!(t, h_atm, Q4278'; rasterize = true)
cb4278 = Colorbar(ga[1, 2], hm4278; label = "photons/m³/s")
ylims!(nothing, 400)
colgap!(ga, 10)
f
