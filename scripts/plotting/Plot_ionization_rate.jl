using AURORA
using MAT
using CairoMakie
CairoMakie.activate!()
using GLMakie
GLMakie.activate!()

set_theme!(fontsize = 20)

## Directory to plot, absolute path
full_path_to_directory = joinpath(REVONTULI_MOUNT,
                                  "mnt/data/etienne/Julia/AURORA.jl/data/Visions2/" *
                                  "InvertedV_480s_fixed-secondaries")

# Read the ionization data
ionization_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")
data = matread(ionization_file)
QO2i = data["QO2i"] .* 1e4 # Q was calculated with Ie in (/cm²) when it sould have been in (/m²), so we fix it here by multiplying with 1e4
QOi = data["QOi"] .* 1e4
QN2i = data["QN2i"] .* 1e4
Qtot = QN2i + QO2i + QOi
h_atm = vec(data["h_atm"])
t = vec(data["t"])


## Plot ionization rate as an animated lineplot
#region
fig = Figure(size=(800, 1200))
# make slider to show and control i_t
sl_time = Slider(fig[2, 1], range = 1:length(t), startvalue = 1)
# make observables
i_t = lift(sl_time.value) do Int; Int; end
time = Observable(string(round(t[i_t[]], digits=3)) * "s")
QO2i_slice = Observable(QO2i[:, i_t[]])
QOi_slice = Observable(QOi[:, i_t[]])
QN2i_slice = Observable(QN2i[:, i_t[]])
Qtot_slice = Observable(QO2i_slice[] + QOi_slice[] + QN2i_slice[])
# plot ionization rate
ax1 = Axis(fig[1, 1], title = time,
           xlabel = "Ionization rate (/m³/s)", ylabel = "altitude (km)",
           yminorticksvisible = false, yminorgridvisible = false, yticks = 100:100:600,
           xminorticksvisible = true, xminorgridvisible = true,
           xminorticks = IntervalsBetween(9),
           xscale = log10, xticks = LogTicks(2:15))
Q_max = maximum(maximum.([QO2i, QOi, QN2i]))
xlims!(ax1, 1e3, Q_max * 10)
ylims!(h_atm[1] / 1e3 - 30, h_atm[end] / 1e3 + 30)
l_QO2i = lines!(QO2i_slice, h_atm / 1e3; linestyle =:dash)
l_QOi = lines!(QOi_slice, h_atm / 1e3; linestyle =:dash)
l_QN2i = lines!(QN2i_slice, h_atm / 1e3; linestyle =:dash)
l_tot = lines!(Qtot_slice, h_atm / 1e3; color = :black)
Legend(fig[1, 2], [l_QO2i, l_QOi, l_QN2i, l_tot], ["O2", "O", "N2", "tot"],
       "Ionization rates")
# make button to run/stop the animation
fig[3, 1] = buttongrid = GridLayout(tellheight = true, tellwidth = false)
run = buttongrid[1, 1] = Button(fig; label = "run")
isrunning = Observable(false)
empty!(run.clicks.listeners)
on(run.clicks) do clicks; isrunning[] = !isrunning[]; end
on(run.clicks) do clicks
    @async while isrunning[]
        isopen(fig.scene) || break # ensures computations stop if closed window
        sl_time.value[] < length(t) ? sl_time.value[] += 1 : sl_time.value[] = 1 # take next time step and loop when t_max is reached
        Makie.set_close_to!(sl_time, sl_time.value[])
        sleep(0.01)
    end
end
# make button to reset animation
reset = buttongrid[1, 2] = Button(fig; label = "reset")
on(reset.clicks) do clicks;
    sl_time.value[] = 1;
    Makie.set_close_to!(sl_time, sl_time.value[])
end
# make function to step in time
function step!(QO2i_slice, QOi_slice, QN2i_slice, Qtot_slice, time, i_t)>
    QO2i_slice[] = QO2i[:, i_t[]]
    QOi_slice[] = QOi[:, i_t[]]
    QN2i_slice[] = QN2i[:, i_t[]]
    Qtot_slice[] = QO2i_slice[] + QOi_slice[] + QN2i_slice[]
    time[] = string(round(t[i_t], digits=3)) * "s"
end
onany(i_t) do i_t
    step!(QO2i_slice, QOi_slice, QN2i_slice, Qtot_slice, time, i_t[])
end
display(GLMakie.Screen(), fig)

## Save animation as a video
video_file = joinpath(full_path_to_directory, "Ionization_rate.mp4")
record(fig, video_file, 1:length(t); framerate = 60, px_per_unit = 2) do i_t
    Makie.set_close_to!(sl_time, i_t)
end
println("Saved $video_file")
#endregion





## Plot ionization rate as heatmap + top incoming flux
#region
# File with Ie_top flux datafilename
# Ietop_file = joinpath(full_path_to_directory, "Ie_incoming_475s.mat")
incoming_files = filter(file -> startswith(file, "Ie_incoming_"), readdir(full_path_to_directory))
if length(incoming_files) > 1
    error("More than one file contains incoming flux. This is not normal")
else
    global Ietop_file = joinpath(full_path_to_directory, incoming_files[1])
end

# Read the Ie_top flux
data = matread(Ietop_file)
Ietop = data["Ie_total"]
θ_Visions_lims = data["theta_lims"]
t_top = data["t_top"]; t_top = [t_top; t_top[end] + diff(t_top)[end]] .- t_top[1]
Ie_file = joinpath(full_path_to_directory, "IeFlickering-01.mat")
Egrid = matopen(Ie_file) do file; read(file, "E"); end # extract Egrid from the simulation
dEgrid = diff(Egrid); dEgrid = [dEgrid; dEgrid[end] + diff(dEgrid)[end]] # calculate dE
# Plot e- precipitation
fig = Figure(size = (850, 950))
angle_cone = [90 180] # angle for the cone of precipitation to plot, 180° is field-aligned down
ax_Ietop = Axis(fig[1, 1], yscale = log10, ylabel = " Energy (eV)",
                title = string(180 - angle_cone[2], " - ", 180 - angle_cone[1], "° DOWN"),
                yminorticksvisible = true, yminorticks = IntervalsBetween(9),
                xticklabelsvisible = false, xminorticksvisible = true,
                xticksmirrored = true, yticksmirrored = true,
                limits = ((0, 1), nothing))
idx_θ = vec(angle_cone[1] .<= abs.(acosd.(mu_avg(θ_Visions_lims))) .<= angle_cone[2])
BeamW = beam_weight([angle_cone[1], angle_cone[2]])
data_heatmap = dropdims(sum(Ietop[idx_θ, :, :]; dims=1); dims=1) ./ BeamW ./ dEgrid' .* Egrid'
hm = Makie.heatmap!(t_top, Egrid, data_heatmap; colorrange = (1e6, maximum(data_heatmap)),
                    colorscale = log10, colormap = :inferno)
custom_formatter(values) = map(v -> rich("10", superscript("$(round(Int64, v))")), values)
Colorbar(fig[1, 2], hm; label = "IeE (eV/m²/s/eV/ster)")
# Plot ionization rate as heatmap
ax_Q = Axis(fig[2, 1], xlabel = "t (s)", ylabel = "altitude (km)", aspect = AxisAspect(1),
          xminorticksvisible = true, yminorticksvisible = true, xticksmirrored = true,
          yticksmirrored = true)
hm = heatmap!(t, h_atm / 1e3, Qtot'; colorrange = (1e6, maximum(Qtot)), colorscale=log10, rasterize = true)
cb = Colorbar(fig[2, 2], hm; label = "Ionization rate (/m³/s)")
Qtot_contour = log10.(Qtot) |> (x -> replace(x, -Inf => 0))
ct = contour!(t, h_atm / 1e3, Qtot_contour'; levels = 6:10, color = :black, labels = true)
colsize!(fig.layout, 1, Aspect(2, 1.0))
rowsize!(fig.layout, 1, Relative(1/5))
linkxaxes!(ax_Ietop, ax_Q)
xlims!(ax_Q, 0, 1)
ylims!(ax_Q, 100, 500)
# display(fig, backend = GLMakie)
display(fig)

## Save the figure
savefile = joinpath(full_path_to_directory, "Ionization_rate_zt.png")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")
savefile = joinpath(full_path_to_directory, "Ionization_rate_zt.svg")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")
savefile = joinpath(full_path_to_directory, "Ionization_rate_zt.pdf")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")
savefile = joinpath(full_path_to_directory, "Ionization_rate_zt.eps")
save(savefile, fig; backend = CairoMakie)
println("Saved $savefile")
#endregion
