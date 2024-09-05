#=
This is a script to produce figures related to the optical aurora task
=#
using AURORA
using MAT
# using CairoMakie
# CairoMakie.activate!()
using GLMakie
GLMakie.activate!()


## Directory to plot, absolute path
# full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/Optical_Aurora_course/input_from_ketchup_3keV-10"
# full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/Optical_Aurora_course/tests/test_constant_onset_20keV_zmin90km_dz15m-from150km"
# full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/Optical_Aurora_course/tests/test_constant_onset_20keV_zmin90km_dz25m-from150km"
# full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/Optical_Aurora_course/tests/test_constant_onset_20keV_zmin90km_dz50m-from150km"
# full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/Optical_Aurora_course/tests/test_constant_onset_20keV_zmin90km_dz75m-from150km"
# full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/Optical_Aurora_course/tests/test_constant_onset_20keV_zmin90km_dz100m-from150km"
# full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/Optical_Aurora_course/tests/test_constant_onset_20keV_zmin90km_dz150m-from150km"
# full_path_to_directory = "data/Optical_Aurora_course/Maxwellian_3kev_with_LET"
# full_path_to_directory = "data/Optical_Aurora_course/Maxwellian_3kev_with_LET_dzchangefrom200km"
full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/Optical_Aurora_course/tests/test_constant_onset_20keV_zmin90km_dz15m-from150km/"
# Read the ionization data
ionization_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")
data = matread(ionization_file)
QO2i = data["QO2i"]
QOi = data["QOi"]
QN2i = data["QN2i"]
Qtot = QN2i + QO2i + QOi
QO1D = data["QO1D"]
h_atm = vec(data["h_atm"]) / 1e3
# t = vec(data["t"])
t = 1



## Plots comparing total ionization profiles for different dz
#region
f = Figure(size = (1000, 1000))
ax = Axis(f[1, 1]; xscale = log10, limits = ((1e6, 1e14), (90, 200)),
          xminorgridvisible = true, xminorticks = IntervalsBetween(9),
          xticks = LogTicks(0:17), title = "")
l_tot_75m = lines!(vec(Qtot) .+ 1e-10, h_atm)
l_N2_75m = lines!(vec(QN2i) .+ 1e-10, h_atm)
l_O2_75m = lines!(vec(QO2i) .+ 1e-10, h_atm)
l_O_75m = lines!(vec(QOi) .+ 1e-10, h_atm)
# Legend(f[1, 2], [l_tot, l_N2, l_O2, l_O], ["tot", "N2", "O2", "O"],
#        "Ionization rates")
display(GLMakie.Screen(), f)

##
l_tot_100m = lines!(vec(Qtot) .+ 1e-10, h_atm; color=Cycled(1), linestyle=:dashdot)
l_N2_100m = lines!(vec(QN2i) .+ 1e-10, h_atm; color=Cycled(2), linestyle=:dashdot)
l_O2_100m = lines!(vec(QO2i) .+ 1e-10, h_atm; color=Cycled(3), linestyle=:dashdot)
l_O_100m = lines!(vec(QOi) .+ 1e-10, h_atm; color=Cycled(4), linestyle=:dashdot)
##
l_tot_150m = lines!(vec(Qtot) .+ 1e-10, h_atm; color=Cycled(1), linestyle=:dot)
l_N2_150m = lines!(vec(QN2i) .+ 1e-10, h_atm; color=Cycled(2), linestyle=:dot)
l_O2_150m = lines!(vec(QO2i) .+ 1e-10, h_atm; color=Cycled(3), linestyle=:dot)
l_O_150m = lines!(vec(QOi) .+ 1e-10, h_atm; color=Cycled(4), linestyle=:dot)
##
# Leg = Legend(f[1, 2], [[l_tot_75m, l_N2_75m, l_O2_75m, l_O_75m], [l_tot_100m, l_N2_100m, l_O2_100m, l_O_100m], [l_tot_150m, l_N2_150m, l_O2_150m, l_O_150m]],
#                 [["tot", "N2", "O2", "O"], ["tot", "N2", "O2", "O"], ["tot", "N2", "O2", "O"]],
#                 [["dz = 150m"], ["dz = 75m"], ["dz = 25m"]])
Leg = Legend(f[1, 2], [[l_tot_75m, l_N2_75m, l_O2_75m, l_O_75m], [l_tot_100m, l_N2_100m, l_O2_100m, l_O_100m]],
                [["tot", "N2", "O2", "O"], ["tot", "N2", "O2", "O"]],
                [["dz = 25m"], ["dz = 15m"]])
#endregion



## Plot the flux as a heatmap in height and energy
#region
data = matread("data/Optical_Aurora_course/input_from_ketchup_3keV-10/IeFlickering-01.mat")
# data = matread("data/Optical_Aurora_course/test_constant_onset_40keV/IeFlickering-01.mat")
# data = matread("/run/user/1000/gvfs/sftp:host=revontuli.uit.no/mnt/data/etienne/MATLAB/AURORA/MI_coupling/data_archives/thesis_with_bug/conversion_1.27e7_120s-4/MIC-18streams-0.35s-4/IeFlickering-01.mat")
μ_lims = data["mu_lims"]
h_atm = data["h_atm"]
E = data["E"]
dE = diff(E); dE = [dE; dE[end]]
Ie = data["Ie_ztE"]


## Plot the Ie over altitude and energy, with different panels for different pitch-angles
using GLMakie
θ_plot =   [(170, 180)  (150, 170) (120, 150)  (90, 120);  # DOWN
            (0, 10)  (10, 30) (30, 60) (60, 90)]  # UP

f = Figure(size = (1800, 1000))
data_heatmap = [zeros(length(h_atm), size(Ie, 2), size(Ie, 3)) for _ in θ_plot]
for i_r in axes(θ_plot, 1)
    for i_c in axes(θ_plot, 2)
        idx_θ = vec(θ_plot[i_r, i_c][1] .<= abs.(acosd.(mu_avg(180:-10:0))) .<= θ_plot[i_r, i_c][2])
        # println(idx_θ)
        BeamW = beam_weight([θ_plot[i_r, i_c][1], θ_plot[i_r, i_c][2]])
        for stream_to_sum in findall(idx_θ)
            idx_z = (1:length(h_atm)) .+ (stream_to_sum - 1) * length(h_atm)
            data_heatmap[i_r, i_c] += Ie[idx_z, :, :] ./ BeamW ./ reshape(dE, (1, 1, :)) #.* reshape(E, (1, 1, :))
        end
        ax = Axis(f[i_r, i_c]; xscale = log10)
        hm = heatmap!(E, h_atm / 1e3, log10.(data_heatmap[i_r, i_c][:, 1, :])';
                      colorrange = (0, 2),
                    #   colorrange = (5, 10),
                      )
        # l = lines!(E, data_heatmap[i_r, i_c][end - 1, :, :][1, :])
        xlims!(2e0, 5e3)
        if i_c == size(θ_plot, 2)
            cb = Colorbar(f[i_r, i_c + 1], hm; label = "log10(#e-/m²/s/eV/ster)")
        end
        if i_r < size(θ_plot, 1)
            ax.title = string(180 - θ_plot[i_r, i_c][2], " - ", 180 - θ_plot[i_r, i_c][1], "° DOWN")
            ax.xlabel = ""
            ax.xticklabelsvisible = false
        else
            ax.title = string(θ_plot[i_r, i_c][1], " - ", θ_plot[i_r, i_c][2], "° UP")
            ax.xlabel = "Energy (eV)"
        end
        if i_c > 1
            ax.yticklabelsvisible = false
        else
            ax.ylabel = "height (km)"
        end
    end
end


display(GLMakie.Screen(), f)


# extremas = filter.(!in(-Inf), data_heatmap) |> (x -> map(extrema, x)) # filtering for the -Inf
# x = maximum(maximum.(extremas))
#endregion



## Plots comparing optical emission profiles for different simulation parameters
#=
Some notes
- The factor 5.627e-3 has to do with taking quenching properly into account (ask Björn)
- The factor 2.3 is because I accidentaly used two different total energy flux for the 1keV
and 3keV precipitation spectra when it should have been exactly the same total energy flux.
=#
#region
# 1st simulation to compare
full_path_to_directory = "data/Optical_Aurora_course/Maxwellian_3kev_with_LET_dzchangefrom200km"
# Read the ionization data
ionization_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")
data = matread(ionization_file)
QO1D = data["Q4278"]
# QO1D = data["QO1D"] .* 5.627e-3
ionization_file = joinpath(full_path_to_directory, "Qzt_all_L_old.mat")
data = matread(ionization_file)
QO1D_old = data["QO1D"]
h_atm = vec(data["h_atm"]) / 1e3
t = 1
# Plot
f = Figure(size = (1000, 1000))
ax = Axis(f[1, 1]; xscale = log10, limits = ((1e-5, 1e4), (90, 400)),
          xminorgridvisible = true, xminorticks = IntervalsBetween(9),
          xticks = LogTicks(-5:17), title = "4278 Å", xlabel = "Emission rate (#/m³/s)", ylabel = "height (km)")
l_O1D_3keV_M = lines!(vec(QO1D) .+ 1e-10, h_atm)
# l_O1D_3keV_M_noquench = lines!(vec(QO1D_old) .+ 1e-10, h_atm, linewidth=1, color=Cycled(1))
display(GLMakie.Screen(), f)

# 2nd simulation to compare
full_path_to_directory = "data/Optical_Aurora_course/Maxwellian_1kev_with_LET_dzchangefrom200km"
# Read the ionization data
ionization_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")
data = matread(ionization_file)
QO1D = data["Q4278"] * 2.3
# QO1D = data["QO1D"] .* 5.627e-3 * 2.3
ionization_file = joinpath(full_path_to_directory, "Qzt_all_L_old.mat")
data = matread(ionization_file)
QO1D_old = data["QO1D"] * 2.3
h_atm = vec(data["h_atm"]) / 1e3
t = 1
# Plot
l_O1D_1keV_M = lines!(vec(QO1D) .+ 1e-10, h_atm)
# l_O1D_1keV_M_noquench = lines!(vec(QO1D_old) .+ 1e-10, h_atm, linewidth=1, color=Cycled(2))

# 3rd simulation to compare
full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/Optical_Aurora_course/input_from_ketchup_3keV-10"
# Read the ionization data
ionization_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")
data = matread(ionization_file)
QO1D = data["Q4278"]
# QO1D = data["QO1D"] .* 5.627e-3
ionization_file = joinpath(full_path_to_directory, "Qzt_all_L_old.mat")
data = matread(ionization_file)
QO1D_old = data["QO1D"]
h_atm = vec(data["h_atm"]) / 1e3
t = 1
# Plot
l_O1D_3keV_G = lines!(vec(QO1D) .+ 1e-10, h_atm; color=Cycled(1), linestyle=:dashdot)
# l_O1D_3keV_G_noquench = lines!(vec(QO1D_old) .+ 1e-10, h_atm; color=Cycled(1), linestyle=:dashdot, linewidth=1)

# 4th simulation to compare
full_path_to_directory = "/home/etienne/Documents/Julia/AURORA.jl/data/Optical_Aurora_course/input_from_ketchup_1keV-9"
# Read the ionization data
ionization_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")
data = matread(ionization_file)
QO1D = data["Q4278"] * 2.3
# QO1D = data["QO1D"] .* 5.627e-3 * 2.3
ionization_file = joinpath(full_path_to_directory, "Qzt_all_L_old.mat")
data = matread(ionization_file)
QO1D_old = data["QO1D"] * 2.3
h_atm = vec(data["h_atm"]) / 1e3
t = 1
l_O1D_1keV_G = lines!(vec(QO1D) .+ 1e-10, h_atm; color=Cycled(2), linestyle=:dashdot)
# l_O1D_1keV_G_noquench = lines!(vec(QO1D_old) .+ 1e-10, h_atm; color=Cycled(2), linestyle=:dashdot, linewidth=1)

# Leg = Legend(f[1, 2], [[l_O1D_3keV_M, l_O1D_3keV_M_noquench, l_O1D_3keV_G, l_O1D_3keV_G_noquench], [l_O1D_1keV_M, l_O1D_1keV_M_noquench, l_O1D_1keV_G, l_O1D_1keV_G_noquench]],
#                 [["Maxwellian emission", "Maxwellian production", "Gaussian emission",  "Gaussian production"], ["Maxwellian emission", "Maxwellian production", "Gaussian emission", "Gaussian production"]],
#                 [["3 keV"], ["1 keV"]])
Leg = Legend(f[1, 2], [[l_O1D_3keV_M, l_O1D_3keV_G], [l_O1D_1keV_M, l_O1D_1keV_G]],
                [["Maxwellian ", "Gaussian "], ["Maxwellian ", "Gaussian "]],
                [["3 keV"], ["1 keV"]])
#endregion




## Checking dz inhomogeneity by plotting profiles of Ie[:, :, iE]
dzchange = 130
file_to_read = "/home/etienne/Documents/Julia/AURORA.jl/data/Optical_Aurora_course/Maxwellian_3kev_with_LET_dzchangefrom130km/IeFlickering-01.mat"
data = matread(file_to_read)
μ_lims = data["mu_lims"]
h_atm = data["h_atm"]
E = data["E"]
dE = diff(E); dE = [dE; dE[end]]
Ie = data["Ie_ztE"]
# Resize Ie from [n_μ x n_z, 1, n_E] to [n_μ, n_z, n_E]
Ie_restructured = zeros(length(μ_lims) - 1, length(h_atm), length(E))
for i in 1:(length(μ_lims) - 1)
    Ie_restructured[i, :, :] .= Ie[(1:length(h_atm)) .+ (i - 1) * length(h_atm), 1, :]
end
# Plot
f = Figure(size = (1800, 1000))
ax = Axis(f[1, 1]; xscale = log10)
iE = findmin(abs.(E .- 10000))[2]
l = lines!(ax, Ie_restructured[1, :, iE], h_atm / 1e3)
iE = findmin(abs.(E .- 5000))[2]
l = lines!(ax, Ie_restructured[1, :, iE], h_atm / 1e3)
iE = findmin(abs.(E .- 3000))[2]
l = lines!(ax, Ie_restructured[1, :, iE], h_atm / 1e3)
iE = findmin(abs.(E .- 1000))[2]
l = lines!(ax, Ie_restructured[1, :, iE], h_atm / 1e3)
iE = findmin(abs.(E .- 500))[2]
l = lines!(ax, Ie_restructured[1, :, iE], h_atm / 1e3)
xlims!(ax, 1e-4, 1e1)
Legend(f[1, 2], ax.scene.plots, ["10 keV", "5 keV", "3 keV", "1 keV", "500 eV"])
hlines!(dzchange; color = :black, linestyle = :dash)
text!(1e-4, dzchange, text = "$dzchange km", align = (:left, :bottom), offset = (10, 0))
display(GLMakie.Screen(), f)




## Checking dz inhomogeneity by plotting profiles of Ie[:, :, iE], comparing simulations
# Load 1st simulation
file_to_read = "/home/etienne/Documents/Julia/AURORA.jl/data/Optical_Aurora_course/tests/test_constant_onset_20keV_zmin90km_nohalfsteps_dz150m-from150km/IeFlickering-01.mat"
data = matread(file_to_read)
μ_lims = data["mu_lims"]
h_atm = data["h_atm"]
E = data["E"]
dE = diff(E); dE = [dE; dE[end]]
Ie = data["Ie_ztE"]
# Resize Ie from [n_μ x n_z, 1, n_E] to [n_μ, n_z, n_E]
Ie_restructured = zeros(length(μ_lims) - 1, length(h_atm), length(E))
for i in 1:(length(μ_lims) - 1)
    Ie_restructured[i, :, :] .= Ie[(1:length(h_atm)) .+ (i - 1) * length(h_atm), 1, :]
end
# Plot
f = Figure(size = (1800, 1000))
ax = Axis(f[1, 1]; xscale = log10)
iE = findmin(abs.(E .- 20000))[2]
l1 = lines!(ax, Ie_restructured[1, :, iE], h_atm / 1e3)
iE = findmin(abs.(E .- 10000))[2]
l2 = lines!(ax, Ie_restructured[1, :, iE], h_atm / 1e3)
iE = findmin(abs.(E .- 5000))[2]
l3 = lines!(ax, Ie_restructured[1, :, iE], h_atm / 1e3)
iE = findmin(abs.(E .- 1000))[2]
l4 = lines!(ax, Ie_restructured[1, :, iE], h_atm / 1e3)
xlims!(ax, 1e-4, 1e1)

# Load 2nd simulation
file_to_read = "/home/etienne/Documents/Julia/AURORA.jl/data/Optical_Aurora_course/tests/test_constant_onset_20keV_zmin90km_nohalfsteps_dz25m-from150km/IeFlickering-01.mat"
data = matread(file_to_read)
μ_lims = data["mu_lims"]
h_atm = data["h_atm"]
E = data["E"]
dE = diff(E); dE = [dE; dE[end]]
Ie = data["Ie_ztE"]
# Resize Ie from [n_μ x n_z, 1, n_E] to [n_μ, n_z, n_E]
Ie_restructured = zeros(length(μ_lims) - 1, length(h_atm), length(E))
for i in 1:(length(μ_lims) - 1)
    Ie_restructured[i, :, :] .= Ie[(1:length(h_atm)) .+ (i - 1) * length(h_atm), 1, :]
end
# Plot
iE = findmin(abs.(E .- 20000))[2]
lines!(ax, Ie_restructured[1, :, iE], h_atm / 1e3; color = l1.color, linestyle = :dashdot)
iE = findmin(abs.(E .- 10000))[2]
lines!(ax, Ie_restructured[1, :, iE], h_atm / 1e3; color = l2.color, linestyle = :dashdot)
iE = findmin(abs.(E .- 5000))[2]
lines!(ax, Ie_restructured[1, :, iE], h_atm / 1e3; color = l3.color, linestyle = :dashdot)
iE = findmin(abs.(E .- 1000))[2]
lines!(ax, Ie_restructured[1, :, iE], h_atm / 1e3; color = l4.color, linestyle = :dashdot)
Legend(f[1, 2], ax.scene.plots, ["20 keV - dz 150m", "10 keV - dz 150m", "5 keV - dz 150m", "1 keV - dz 150m",
                                 "20 keV - dz 25m", "10 keV - dz 25m", "5 keV - dz 25m", "1 keV - dz 25m"])
display(GLMakie.Screen(), f)

xlims!(1e4, 1e11)
ylims!(80, 200)
