## SHOULD WE TAKE INTO ACCOUNT THERMAL E- ???

using AURORA
using MAT
using GLMakie
# using CairoMakie
# CairoMakie.activate!()
# using WGLMakie
# WGLMakie.activate!()


directory_to_plot = "Visions2/Alfven_536s_correct_msis_and_scattering"



## Read the data
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
i_t = Observable(1)
E_limit = 5
n_e_superthermal = Observable(dropdims(sum(n_e[:, i_t[], :], dims=2), dims=2))
n_e_over100eV = Observable(dropdims(sum(n_e[:, i_t[], E .> E_limit], dims=2), dims=2))
n_e_under100eV = Observable(dropdims(sum(n_e[:, i_t[], E .< E_limit], dims=2), dims=2))
ne_total = Observable(ne_background .+ n_e_superthermal[])
time = Observable(string(round(t[i_t[]], digits=3)) * "s")

fig = Figure()
custom_formatter(values) = map(v -> "10" * Makie.UnicodeFun.to_superscript(round(Int64, v)), values)
ax1 = Axis(fig[1, 1], title = time, xlabel = "nₑ (m⁻³)", ylabel = "altitude (km)",
    yminorticksvisible = false, yminorgridvisible = false, yticks = 100:100:600,
    xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(9),
    xscale = log10, xticks = LogTicks(4:8)
)
n_e_max = maximum(sum(n_e, dims=3))
xlims!(ax1, 1e3, n_e_max * 10)
ylims!(50, h_atm[end] / 1e3 + 50)
l_superthermal = lines!(n_e_superthermal, h_atm / 1e3)
l_over100ev = lines!(n_e_over100eV, h_atm / 1e3)
l_under100ev = lines!(n_e_under100eV, h_atm / 1e3)
Legend(fig[1, 2], [l_superthermal, l_over100ev, l_under100ev], ["total superthermal", ">  $E_limit eV", "< $E_limit eV"])

ax2 = Axis(fig[1, 4], title = time, xlabel = "nₑ (m⁻³)", ylabel = "altitude (km)",
    yminorticksvisible = false, yminorgridvisible = false, yticks = 100:100:600,
    xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(9),
    xscale = log10,
)
ne_max = maximum(ne_background) + n_e_max
xlims!(ax2, 1e9, ne_max * 10)
ylims!(50, h_atm[end] / 1e3 + 50)
l_background = lines!(ne_background, h_atm / 1e3)
l_tot = lines!(ne_total, h_atm / 1e3; linestyle=:dash)

Legend(fig[1, 5], [l_background, l_tot], ["background", "background + superthermal"])
display(fig)

## try to make that interactive
function step!(n_e_superthermal, n_e_over100eV, n_e_under100eV, ne_total, time, i_t)
    n_e_superthermal[] = dropdims(sum(n_e[:, i_t, :], dims=2), dims=2)
    n_e_over100eV[] = dropdims(sum(n_e[:, i_t, E .> E_limit], dims=2), dims=2)
    n_e_under100eV[] = dropdims(sum(n_e[:, i_t, E .< E_limit], dims=2), dims=2)
    ne_total[] = ne_background .+ n_e_superthermal[]

    time[] = string(round(t[i_t], digits=3)) * "s"
end

# make button to run/stop the animation
fig[1, 3] = buttongrid = GridLayout(tellheight = false, tellwidth = false)
run = buttongrid[1, 1] = Button(fig; label = "run")
isrunning = Observable(false)
empty!(run.clicks.listeners)
on(run.clicks) do clicks; isrunning[] = !isrunning[]; end
on(run.clicks) do clicks
    @async while isrunning[]
        isopen(fig.scene) || break # ensures computations stop if closed window
        i_t[] < length(t) ? i_t[] += 1 : i_t[] = 1 # take next time step and loop when t_max is reached
        step!(n_e_superthermal, n_e_over100eV, n_e_under100eV, ne_total, time, i_t[])
        sleep(0.01)
    end
end
reset = buttongrid[2, 1] = Button(fig; label = "reset")
on(reset.clicks) do clicks; i_t[] = 1; end
