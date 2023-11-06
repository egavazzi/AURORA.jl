using AURORA
using MAT
using CairoMakie
CairoMakie.activate!()
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

## Plot e- density
# i_t = 100
# n_e_to_plot = dropdims(sum(n_e[:, i_t, :], dims=2), dims=2)

n_e_to_plot = Observable(dropdims(sum(n_e[:, 1, :], dims=2), dims=2))

f = Figure()
ax = Axis(f[1, 1])
l = lines!(n_e_to_plot, h_atm / 1e3)
f

for i_t in eachindex(t)
    n_e_to_plot[] = dropdims(sum(n_e[:, i_t, :], dims=2), dims=2)
    sleep(1/30)
    display(f)
end
