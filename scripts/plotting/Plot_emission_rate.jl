# This script produces plots of altitude-time-variation of volume emission-rates and
# time-variation of normalized column-emission (excitation) intensity plots.

using AURORA
using MAT
using CairoMakie
using GLMakie
GLMakie.activate!()

## Directory to plot, absolute path
full_path_to_directory = joinpath(REVONTULI_MOUNT,
                                  "mnt/data/etienne/Julia/AURORA.jl/data/Visions2/" *
                                  "InvertedV_480s_fixed-secondaries")


## Load the Q data (volume emission-rates)
Q_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")
data = matread(Q_file)
Q4278 = data["Q4278"] * 1e4
Q6730 = data["Q6730"] * 1e4
Q7774 = data["Q7774"] * 1e4
Q8446 = data["Q8446"] * 1e4
QO1D = data["QO1D"] * 1e4
QO1S = data["QO1S"] * 1e4
h_atm = vec(data["h_atm"]) ./ 1e3 # convert to km
t = vec(data["t"])

## Plot Q data
fig = with_theme(
    Theme(
        Axis = (
            xticksmirrored = true, yticksmirrored = false, xminorticksvisible = true,
            yminorticksvisible = true, limits=((0, 1), (100, 400))
            ),
        Heatmap = (rasterize = true,),
        Colorbar = (flip_vertical_label = true, vertical = true),
        ),
        fontsize = 20
    ) do
    fig = Figure(size = (1000, 800))
    ga = fig[1, 1] = GridLayout()
    ax4278 = Axis(ga[1, 1]; title = "4278 Å", xticklabelsvisible = false, ylabel ="altitude (km)")
    hm4278 = heatmap!(t, h_atm, Q4278'; rasterize = true)
    cb4278 = Colorbar(ga[1, 2], hm4278; label = "photons/m³/s")
    colgap!(ga, 10)

    gb = fig[1, 2] = GridLayout()
    ax6730 = Axis(gb[1, 1]; title = "6730 Å", xticklabelsvisible = false, yticklabelsvisible = false)
    hm6730 = heatmap!(t, h_atm, Q6730')
    cb6730 = Colorbar(gb[1, 2], hm6730; label = "photons/m³/s")
    colgap!(gb, 10)

    gc = fig[2, 1] = GridLayout()
    ax7774 = Axis(gc[1, 1]; title = "7774 Å", xlabel = "time(s)", ylabel = "altitude (km)")
    hm7774 = heatmap!(t, h_atm, Q7774')
    cb7774 = Colorbar(gc[1, 2], hm7774; label = "photons/m³/s")
    colgap!(gc, 10)

    gd = fig[2, 2] = GridLayout()
    ax8446 = Axis(gd[1, 1]; title = "8446 Å", yticklabelsvisible = false, xlabel = "time (s)")
    hm8446 = heatmap!(t, h_atm, Q8446')
    cb8446 = Colorbar(gd[1, 2], hm8446; label = "photons/m³/s")
    colgap!(gd, 10)
    return fig
end
display(fig)

## Save Qtz plot
savefile = joinpath(full_path_to_directory, "Qtz.png")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")
savefile = joinpath(full_path_to_directory, "Qtz.svg")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")
savefile = joinpath(full_path_to_directory, "Qtz.pdf")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")
savefile = joinpath(full_path_to_directory, "Qtz.eps")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")
##





## Load the I data (column-integrated emission-rates)
I_file = joinpath(full_path_to_directory, "I_lambda_of_t.mat")
data = matread(I_file)
I4278 = vec(data["I_4278"]) * 1e4 / 1e10
I6730 = vec(data["I_6730"]) * 1e4 / 1e10
I7774 = vec(data["I_7774"]) * 1e4 / 1e10
I7774_O = vec(data["I_7774_O"]) * 1e4 / 1e10
I7774_O2 = vec(data["I_7774_O2"]) * 1e4 / 1e10
I8446 = vec(data["I_8446"]) * 1e4 / 1e10
I8446_O = vec(data["I_8446_O"]) * 1e4 / 1e10
I8446_O2 = vec(data["I_8446_O2"]) * 1e4 / 1e10
IO1D = vec(data["I_O1D"]) * 1e4 / 1e10
IO1S = vec(data["I_O1S"]) * 1e4 / 1e10
t = vec(data["t"])

## Plot the I data
custom_theme = Theme(fontsize = 20, linewidth = 2)
set_theme!(custom_theme)
fig = Figure(size = (1000, 800))
ax = Axis(fig[1, 1]; title = "Intensity", xlabel = "time (s)", ylabel = "R",
          yscale = log10, limits = (t[1], t[end], 1e1, 1e4), xticks = t[1]:0.1:t[end],
          yminorticksvisible = true, yminorgridvisible = true,
          yminorticks = IntervalsBetween(9), yticksmirrored = true)
ylims!(1e1, 1e6)
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
Legend(fig[1, 2], ax; patchsize = [40, 20])
set_theme!()
display(fig)

## Save It plot
savefile = joinpath(full_path_to_directory, "It.png")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")
savefile = joinpath(full_path_to_directory, "It.svg")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")
savefile = joinpath(full_path_to_directory, "It.pdf")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")
savefile = joinpath(full_path_to_directory, "It.eps")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")
##
