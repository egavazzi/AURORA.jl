using AURORA
using MAT
using GLMakie
GLMakie.activate!()


## This is for simple cases of comparison where z-grid, and t-grid are the same
# Path to dta folder 1
# full_path_to_directory1 = joinpath(REVONTULI_MOUNT, "mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_536s_lr_Bz-9_newZ_550km_finer-theta_halfstepsAB")
full_path_to_directory1 = joinpath(REVONTULI_MOUNT, "mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_475s_lr_Bz-9_newZ_550km_finer-theta_halfstepsAB_scaled")
# Path to data folder 2
full_path_to_directory2 = joinpath(REVONTULI_MOUNT, "mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_536s_lr_Bz-9_newZ_550km_finer-theta_halfstepsAB_scaled")


## Read the data
# From data folder 1
ionization_file = joinpath(full_path_to_directory1, "Qzt_all_L.mat")
data = matread(ionization_file)
QO2i = data["QO2i"] .* 1e4 # Q was calculated with Ie in /cm² when it sould have been in /m², so we fix it here by multiplying with 1e4
QOi = data["QOi"] .* 1e4
QN2i = data["QN2i"] .* 1e4
Qtot1 = QN2i + QO2i + QOi
h_atm = vec(data["h_atm"]) / 1e3 # convert to km
t = vec(data["t"])
# From data folder 2
ionization_file = joinpath(full_path_to_directory2, "Qzt_all_L.mat")
data = matread(ionization_file)
QO2i = data["QO2i"] .* 1e4 # Q was calculated with Ie in /cm² when it sould have been in /m², so we fix it here by multiplying with 1e4
QOi = data["QOi"] .* 1e4
QN2i = data["QN2i"] .* 1e4
Qtot2 = QN2i + QO2i + QOi
# h_atm and t are not loaded again (they are the same)

# Resize the ionization matrices dependent on the one with smallest size
if size(Qtot1, 1) > size(Qtot2, 1)
    Qtot1 = Qtot1[1:size(Qtot2, 1), :]
    h_atm = h_atm[1:size(Qtot2, 1)]
elseif size(Qtot1, 1) < size(Qtot2, 1)
    Qtot2 = Qtot2[1:size(Qtot1, 1), :]
    h_atm = h_atm[1:size(Qtot1, 1)]
end

## Plot ionization rates + difference
fig = Figure(size = (1100, 1000))
custom_formatter(values) = map(v -> rich("10", superscript("$(round(Int64, v))")), values)
ga = fig[1, 1] = GridLayout()
ax1 = Axis(ga[1, 1], xlabel = "t (s)", ylabel = "altitude (km)", title = split(full_path_to_directory1, "/")[end])
hm1 = heatmap!(t, h_atm, log10.(Qtot1)'; colorrange = (6, 10))
cb1 = Colorbar(ga[1, 2], hm1; tickformat = custom_formatter, label = "Ionization rate (/m³/s)", tellheight = false)
Qtot_contour1 = log10.(Qtot1) |> (x -> replace(x, -Inf => 0))
ct1 = contour!(t, h_atm, Qtot_contour1'; levels = 6:10, color = :black, labels = true)

gb = fig[2, 1] = GridLayout()
ax2 = Axis(gb[1, 1], xlabel = "t (s)", ylabel = "altitude (km)", title = split(full_path_to_directory2, "/")[end])
hm2 = heatmap!(t, h_atm, log10.(Qtot2)'; colorrange = (6, 10))
cb2 = Colorbar(gb[1, 2], hm2; tickformat = custom_formatter, label = "Ionization rate (/m³/s)")
Qtot_contour2 = log10.(Qtot2) |> (x -> replace(x, -Inf => 0))
ct2 = contour!(t, h_atm, Qtot_contour2'; levels = 6:10, color = :black, labels = true)

gc = fig[1, 2] = GridLayout()
ax3 = Axis(gc[1, 1], xlabel = "t (s)", ylabel = "altitude (km)", title = "ax1 - ax2")
function f(x1, x2)
    y = x1 - x2
    if y > 0
        (y < 1) ? y = 0 : y = log10(y)
    elseif y < 0
        (y > -1) ? y = 0 : y = -log10(abs(y))
    end
    return y
end
Q_compare = f.(Qtot1, Qtot2)
hm3 = heatmap!(t, h_atm, Q_compare'; colormap = :RdBu, colorrange = (-10, 10))
function custom_formatter(values)
    map(v -> begin
            if v > 0
                rich("10", superscript("$(round(Int64, v))"))
            elseif v < 0
                rich("-10", superscript("$(round(Int64, abs(v)))"))
            else # so when v = 0
                "0"
            end
        end, values)
end
# Q_compare = log10.(abs.(Qtot1 .- Qtot2))
# Q_compare = Qtot1 .- Qtot2 .|> (x -> (x < 0) ? 0 : x)
# hm3 = heatmap!(t, h_atm / 1e3, Q_compare'; colormap = :RdBu, colorrange = (-1e9, 1e9))
cb3 = Colorbar(gc[1, 2], hm3; label = "Ionization rate (/m³/s)", tickformat = custom_formatter)

gd = fig[2, 2] = GridLayout()
ax4 = Axis(gd[1, 1], xlabel = "t (s)", ylabel = "altitude (km)", title = "(ax1 - ax2) / ax1")
Q_compare = (Qtot1 .- Qtot2) ./ Qtot1
Q_compare[abs.(Qtot1 .- Qtot2) .< 1e6] .= 0
hm4 = heatmap!(t, h_atm, Q_compare'; colormap = :RdBu, colorrange = (-0.5, 0.5))
cb4 = Colorbar(gd[1, 2], hm4; label = "Relative difference in ionization rate")



display(GLMakie.Screen(), fig)
# display(fig)
