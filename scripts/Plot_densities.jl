## SHOULD WE TAKE INTO ACCOUNT THERMAL E- ???

using AURORA
using MAT
using GLMakie



directory_to_plot = "Visions2/Alfven_536s_correct_msis_and_scattering"


## Density data
# Read the data
full_path_to_directory = pkgdir(AURORA, "data", directory_to_plot)
density_file = joinpath(full_path_to_directory, "superthermal_e_density.mat")
data = matread(density_file)
n_e = data["n_e"]       # [n_z, n_t, n_E]
h_atm = data["h_atm"]   # [n_z]
t = data["t"]           # [n_t]
E = data["E"]           # [n_E]
file_neutral = joinpath(full_path_to_directory, "neutral_atm.mat")
data_neutral = matread(file_neutral)
ne_background = data_neutral["ne"]

## Plot e- density
fig = Figure()
# make slider to show and control i_t
sl_time = Slider(fig[2, 1], range = 1:length(t), startvalue = 1)
# make observables
E_limit = 3
i_t = lift(sl_time.value) do Int; Int; end
time = Observable(string(round(t[i_t[]], digits=3)) * "s")
n_e_superthermal = Observable(dropdims(sum(n_e[:, i_t[], :], dims=2), dims=2))
n_e_over = Observable(dropdims(sum(n_e[:, i_t[], E .> E_limit], dims=2), dims=2))
n_e_under = Observable(dropdims(sum(n_e[:, i_t[], E .< E_limit], dims=2), dims=2))
ratio_over = Observable(n_e_over[] ./ n_e_superthermal[])
ratio_under = Observable(n_e_under[] ./ n_e_superthermal[])
ne_total = Observable(ne_background .+ n_e_superthermal[])
# plot density of superthermal
ax1 = Axis(fig[1, 1], title = time, xlabel = "nₑ (m⁻³)", ylabel = "altitude (km)",
    yminorticksvisible = false, yminorgridvisible = false, yticks = 100:100:600,
    xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(9),
    xscale = log10, xticks = LogTicks(4:8)
)
n_e_max = maximum(sum(n_e, dims=3))
xlims!(ax1, 1e3, n_e_max * 10)
ylims!(50, h_atm[end] / 1e3 + 50)
l_superthermal = lines!(n_e_superthermal, h_atm / 1e3)
l_over = lines!(n_e_over, h_atm / 1e3)
l_under = lines!(n_e_under, h_atm / 1e3)
Legend(fig[1, 2], [l_superthermal, l_over, l_under], ["total superthermal", ">  $E_limit eV", "< $E_limit eV"], "Densities")
# plot density of background + superthermal
ax2 = Axis(fig[1, 3], title = time, xlabel = "nₑ (m⁻³)", ylabel = "altitude (km)",
    yminorticksvisible = false, yminorgridvisible = false, yticks = 100:100:600,
    xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(9),
    xscale = log10,
)
ne_max = maximum(ne_background) + n_e_max
xlims!(ax2, 1e9, ne_max * 10)
ylims!(50, h_atm[end] / 1e3 + 50)
l_background = lines!(ne_background, h_atm / 1e3)
l_tot = lines!(ne_total, h_atm / 1e3; linestyle=:dash)
Legend(fig[1, 4], [l_background, l_tot], ["background", "background + superthermal"], "Densities")
# plot ratio of densities under and over a certain energy
ax3 = Axis(fig[1, 5], title = time, xlabel = "ratio", ylabel = "altitude (km)",
    yminorticksvisible = false, yminorgridvisible = false, yticks = 100:100:600,
    xminorticksvisible = true, xminorgridvisible = true,
)
xlims!(ax3, -0.1, 1.1)
ylims!(50, h_atm[end] / 1e3 + 50)
l_ratio_over = lines!(ratio_over, h_atm / 1e3; color = :green)
l_ratio_under = lines!(ratio_under, h_atm / 1e3; color = :red)
Legend(fig[1, 6], [l_ratio_over, l_ratio_under], ["over $E_limit eV", "under $E_limit eV"], "Density ratio")
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
function step!(n_e_superthermal, n_e_over, n_e_under, ratio_over, ratio_under, ne_total, time, i_t)
    n_e_superthermal[] = dropdims(sum(n_e[:, i_t, :], dims=2), dims=2)
    n_e_over[] = dropdims(sum(n_e[:, i_t, E .> E_limit], dims=2), dims=2)
    n_e_under[] = dropdims(sum(n_e[:, i_t, E .< E_limit], dims=2), dims=2)
    ne_total[] = ne_background .+ n_e_superthermal[]
    ratio_over[] = n_e_over[] ./ n_e_superthermal[]
    ratio_under[] = n_e_under[] ./ n_e_superthermal[]
    time[] = string(round(t[i_t], digits=3)) * "s"
end
onany(i_t) do i_t
    step!(n_e_superthermal, n_e_over, n_e_under, ratio_over, ratio_under, ne_total, time, i_t[])
end
# and display
display(fig)









## Ionization data
# Read the data
full_path_to_directory = pkgdir(AURORA, "data", directory_to_plot)
ionization_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")
data = matread(ionization_file)
QO2i = data["QO2i"] .* 1e4 # the Q were calculated from Ie in m⁻², but we have Ie in cm⁻², so we need to multiply by 1e4
QOi = data["QOi"] .* 1e4
QN2i = data["QN2i"] .* 1e4

## Plot e- density
fig = Figure()
# make slider to show and control i_t
sl_time = Slider(fig[2, 1], range = 1:length(t), startvalue = 1)
# make observables
i_t = lift(sl_time.value) do Int; Int; end
time = Observable(string(round(t[i_t[]], digits=3)) * "s")
QO2i_slice = Observable(QO2i[:, i_t[]])
QOi_slice = Observable(QOi[:, i_t[]])
QN2i_slice = Observable(QN2i[:, i_t[]])
Q_tot = Observable(QO2i_slice[] + QOi_slice[] + QN2i_slice[])
# plot ionization rate
ax1 = Axis(fig[1, 1], title = time, xlabel = "Ionization rate (/m³/s)", ylabel = "altitude (km)",
    yminorticksvisible = false, yminorgridvisible = false, yticks = 100:100:600,
    xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(9),
    xscale = log10, xticks = LogTicks(2:12)
)
Q_max = maximum(maximum.([QO2i, QOi, QN2i]))
xlims!(ax1, 1e3, Q_max * 10)
ylims!(50, h_atm[end] / 1e3 + 50)
l_QO2i = lines!(QO2i_slice, h_atm / 1e3; linestyle =:dash)
l_QOi = lines!(QOi_slice, h_atm / 1e3; linestyle =:dash)
l_QN2i = lines!(QN2i_slice, h_atm / 1e3; linestyle =:dash)
l_tot = lines!(Q_tot, h_atm / 1e3; color = :black)
Legend(fig[1, 2], [l_QO2i, l_QOi, l_QN2i, l_tot], ["O2", "O", "N2", "tot"], "Ionization rates")
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
function step!(QO2i_slice, QOi_slice, QN2i_slice, Q_tot, time, i_t)
    QO2i_slice[] = QO2i[:, i_t[]]
    QOi_slice[] = QOi[:, i_t[]]
    QN2i_slice[] = QN2i[:, i_t[]]
    Q_tot[] = QO2i_slice[] + QOi_slice[] + QN2i_slice[]
    time[] = string(round(t[i_t], digits=3)) * "s"
end
onany(i_t) do i_t
    step!(QO2i_slice, QOi_slice, QN2i_slice, Q_tot, time, i_t[])
end
# and display
display(fig)
