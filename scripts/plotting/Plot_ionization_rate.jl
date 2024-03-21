using AURORA
using MAT
using GLMakie
GLMakie.activate!()




# full_path_to_directory = joinpath(REVONTULI_MOUNT, "mnt/data/etienne/Julia/AURORA.jl/data/steady_state/test_AB-half-steps")
full_path_to_directory = joinpath(REVONTULI_MOUNT, "mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_536s")
# pkgdir(AURORA, "data", "Visions2/Alfven_536s_correct_msis_and_scattering")



## Plot Ionization rate
# Read the data
ionization_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")
data = matread(ionization_file)
QO2i = data["QO2i"] .* 1e4 # the Q were calculated from Ie in m⁻², but we have Ie in cm⁻², so we need to multiply by 1e4
QOi = data["QOi"] .* 1e4
QN2i = data["QN2i"] .* 1e4
h_atm = vec(data["h_atm"])
t = data["t"]


## Plot ionization rates
fig = Figure(size=(800, 1200))
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
    xscale = log10, xticks = LogTicks(2:15)
)
Q_max = maximum(maximum.([QO2i, QOi, QN2i]))
xlims!(ax1, 1e3, Q_max * 10)
ylims!(h_atm[1] / 1e3 - 30, h_atm[end] / 1e3 + 30)
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
function step!(QO2i_slice, QOi_slice, QN2i_slice, Q_tot, time, i_t)>
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





##
savefile = joinpath(full_path_to_directory, "Qz.png")
save(savefile, fig)




##
video_file = joinpath(full_path_to_directory, "Ionization_rate.mp4")
@time record(fig, video_file, 1:length(t); framerate = 60) do i_t
    Makie.set_close_to!(sl_time, i_t)
end
