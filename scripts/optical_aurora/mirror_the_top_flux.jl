using AURORA
using MAT
using GLMakie
GLMakie.activate!()

## Directory to plot, absolute path
full_path_to_directory = "data/Optical_Aurora_course/input_from_ketchup_5keV-10"

# Load the simulation data
data = matread(joinpath(full_path_to_directory, "IeFlickering-01.mat"))
Ie = data["Ie_ztE"]
E = data["E"]
dE = diff(E); dE = diff(E); dE = [dE; dE[end]];
h_atm = data["h_atm"]
μ_lims = data["mu_lims"]
BW = beam_weight(acosd.(μ_lims))
# Restructure Ie from [n_μ x n_z, 1, n_E] to [n_μ, n_z, n_E]
Ie_restructured = zeros(length(μ_lims) - 1, length(h_atm), length(E))
for i_μ in 1:(length(μ_lims) - 1)
    idx_z = (1:length(h_atm)) .+ (i_μ - 1) * length(h_atm)
    Ie_restructured[i_μ, :, :] = Ie[idx_z, 1, :]
end


## Load the incoming data
# data = matread(joinpath(full_path_to_directory, "Ie_incoming.mat"))
data = matread("data/Optical_Aurora_course/input_from_ketchup_1keV-1/Ie_incoming.mat")
Ietop = data["Ie_total"]


## Plot the flux coming DOWN at the top
f = Figure(size = (1600, 1000))
for i_μ in axes(Ietop, 2)
    if i_μ <= Int(size(Ietop, 2) / 2)
        ax = Axis(f[1, i_μ]; xscale = log10)
        l = lines!(E, Ietop[i_μ][1:length(E)] ./ dE ./ BW[i_μ])
    else
        ax = Axis(f[2, 19 - i_μ]; xscale = log10)
        l = lines!(E, Ietop[i_μ][1:length(E)] ./ dE ./ BW[i_μ])
    end
end
display(GLMakie.Screen(), f)


## Plot the flux coming UP at the top
f = Figure(size = (1600, 1000))
for i_μ in axes(Ie_restructured, 1)
    if i_μ <= Int(size(Ie_restructured, 1) / 2)
        ax = Axis(f[1, i_μ]; xscale = log10)
        l = lines!(E, Ie_restructured[i_μ, end, :] ./ dE ./ BW[i_μ])
    else
        ax = Axis(f[2, 19 - i_μ]; xscale = log10)
        l = lines!(E, Ie_restructured[i_μ, end, :] ./ dE ./ BW[i_μ])
    end
end
display(GLMakie.Screen(), f)





## Reflect the e-flux
#=
The Ie_restructured[10:18, :, :] correspond to the e-flux coming up.
The Ietop[1:9] correspond to the e- flux coming down.

The double layer in the ketchup simulation had a potential drop of 3keV. So we will reflect
all e-flux with energy lower than that.
=#

E_reflection = 1000
iE_reflection = findmin(abs.(E .- E_reflection))[2]

Ietop_reflected = deepcopy(Ietop)
for i_μ in 1:9
    Ietop_reflected[i_μ][1:length(E), :] .+= Ie_restructured[19 - i_μ, end, :]
end

## Plot the reflected flux
f = Figure(size = (1600, 1000))
for i_μ in axes(Ietop, 2)
    if i_μ <= Int(size(Ietop_reflected, 2) / 2)
        ax = Axis(f[1, i_μ]; xscale = log10)
        # l = lines!(E, Ietop[i_μ][1:length(E)] ./ dE ./ BW[i_μ])
        l = lines!(E, Ie_restructured[i_μ, end, :] ./ dE ./ BW[i_μ])
        l = lines!(E, Ietop_reflected[i_μ][1:length(E)] ./ dE ./ BW[i_μ])
    else
        ax = Axis(f[2, 19 - i_μ]; xscale = log10)
        # l = lines!(E, Ietop[i_μ][1:length(E)] ./ dE ./ BW[i_μ])
        l = lines!(E, Ie_restructured[i_μ, end, :] ./ dE ./ BW[i_μ])
        l = lines!(E, Ietop_reflected[i_μ][1:length(E)] ./ dE ./ BW[i_μ])
    end
end
display(GLMakie.Screen(), f)



## Save the original + reflected incoming Ie
filename = joinpath(full_path_to_directory, "Ie_reflected.mat")
matopen(filename, "w") do f
    write(f, "Ie_total", Ietop_reflected)
end









## Compare the incoming flux for the different simulations
full_path_to_directory = [
    "data/Optical_Aurora_course/input_from_ketchup_3keV-1/",
    "data/Optical_Aurora_course/input_from_ketchup_3keV-2/",
    # "data/Optical_Aurora_course/input_from_ketchup_3keV-3/",
    "data/Optical_Aurora_course/input_from_ketchup_3keV-4/",
    # "data/Optical_Aurora_course/input_from_ketchup_3keV-5/",
    "data/Optical_Aurora_course/input_from_ketchup_3keV-6/",
    # "data/Optical_Aurora_course/input_from_ketchup_3keV-7/",
    # "data/Optical_Aurora_course/input_from_ketchup_3keV-8/",
    "data/Optical_Aurora_course/input_from_ketchup_3keV-9/",
    "data/Optical_Aurora_course/input_from_ketchup_3keV-10/",
    ]

# Initialize

set_theme!(Theme(fontsize = 20, linewidth=3))
f = Figure(size = (1000, 750))
ax = []
for i in 1:1
    push!(ax, Axis(f[1, i]; xscale = log10, yscale = log10, xminorgridvisible = true,
               yminorgridvisible = true, xminorticksvisible = true, yminorticksvisible = true, xminorticks = IntervalsBetween(9),
               yminorticks = IntervalsBetween(9), xlabel = "Energy (eV)", ylabel = "Ie (#e-/m²/s/eV/ster)"))
    ylims!(ax[i], 1e-1, 1e2)
    ax[i].title = "$(Int(i * 10 - 10)) - $(Int(i * 10))° down"
    ax[i].xticklabelsvisible = true
    if i >= 2
        ax[i].yticklabelsvisible = false
    end
end
# for i in 18:18
#     push!(ax, Axis(f[2, 19 - i]; xscale = log10, yscale = log10))
#     # ylims!(ax[i], 1e-1, 1e2)
#     # ylims!(ax[i], 1e-1, 1e2)
#     # ax[i].title = "$(Int(180 - i * 10)) - $(Int(190 - i * 10))° up"
#     # if i <= 18
#     #     ax[i].yticklabelsvisible = false
#     # end
#     ylims!(ax[2], 1e-1, 1e2)
#     ax[2].title = "$(Int(180 - i * 10)) - $(Int(190 - i * 10))° up"
#     if i <= 18
#         ax[2].yticklabelsvisible = false
#     end
# end

# Loop
l = []
for (idx_simulation, name_simulation) in enumerate(full_path_to_directory)
    data = matread(joinpath(name_simulation, "IeFlickering-01.mat"))
    Ie = data["Ie_ztE"]
    h_atm = data["h_atm"]
    if idx_simulation <= 1
        E = data["E"]
        dE = diff(E); dE = diff(E); dE = [dE; dE[end]];
        μ_lims = data["mu_lims"]
        BW = beam_weight(acosd.(μ_lims))
    end

    Ie_restructured = zeros(length(μ_lims) - 1, length(h_atm), length(E))
    for i_μ in 1:(length(μ_lims) - 1)
        idx_z = (1:length(h_atm)) .+ (i_μ - 1) * length(h_atm)
        Ie_restructured[i_μ, :, :] = Ie[idx_z, 1, :]
    end
    # Then plot
    for i_μ in 1:1#axes(Ie_restructured, 1)
        # l = lines!(ax[i_μ], E, Ie_restructured[i_μ, end, :] ./ dE ./ BW[i_μ])
        push!(l, lines!(ax[i_μ], E, Ie_restructured[i_μ, end, :] ./ dE ./ BW[i_μ]))
    end
end
l[end].linestyle = :dashdot
Legend(f[1, 2], [l...], ["simulation 1", "simulation 2", "simulation 4", "simulation 6", "simulation 9", "simulation 10"])
display(GLMakie.Screen(), f)

set_theme!()
