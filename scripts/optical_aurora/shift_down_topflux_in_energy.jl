#=
This script can be used to shift down in energy the precipitating electron flux spectra
produced by ketchup (Vlasov simulation of potential drop along a magnetic field line).
Let's say you have a Maxwellian spectra with a characteristic energy around 5keV. With this
script you can shift it down to 3keV while keeping the shape of the spectra.
=#
using AURORA
using MAT
using GLMakie
GLMakie.activate!()

## Load the incoming data
data = matread("data/Optical_Aurora_course/ketchup_inputs_non-isotropic/input_from_ketchup_5keV-1/Ie_incoming.mat")
Ietop = data["Ie_total"]

E, dE = make_energy_grid(6000)
lines(E, dE)

# check Ietop
f, ax, l = lines(E, Ietop[1][1:length(E), 1])

idx_old = findmin(abs.(E .- 2400))[2] # 294
idx_new = findmin(abs.(E .- 2400))[2] # 122
diff_idx = idx_old - idx_new # 172
nE = length(Ietop[1]) # 701
width_idx = nE - idx_old # 407
idx_new + width_idx # 529

Ietop_downshifted = deepcopy(Ietop)
for i_μ in 1:9
    Ietop_downshifted[i_μ][idx_new:(idx_new + width_idx), 1] .= Ietop[i_μ][idx_old:701, 1]
    Ietop_downshifted[i_μ][(idx_new + width_idx + 1):end, 1] .= 0
end


# lines!(E .+ 2000, Ietop_downshifted[1][1:length(E), 1])
# f
##
BW = beam_weight(180:-10:0)
f, ax, l = lines(E, Ietop[1][1:length(E), 1] ./ BW[1])
lines!(E, Ietop_downshifted[1][1:length(E), 1] ./ BW[1] ./ 1.5)
lines!(E, Ietop_downshifted[2][1:length(E), 1] ./ BW[2])
lines!(E, Ietop_downshifted[3][1:length(E), 1] ./ BW[3])
f

##
filename =  joinpath("data/Optical_Aurora_course/input_from_ketchup_1keV-1" ,"Ie_incoming.mat")
filename = rename_if_exists(filename)
matopen(filename, "w") do f
    write(f, "Ie_total", Ietop_downshifted)
end
