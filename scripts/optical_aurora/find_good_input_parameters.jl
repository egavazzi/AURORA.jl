using AURORA
using MAT

altitude_max = 500;         # (km) top altitude of the ionosphere
θ_lims = 180:-10:0;         # (°) angle-limits for the electron beams
E_max = 20000;               # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;     # (°) angle between the B-field line and the zenith

msis_file = find_nrlmsis_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );
iri_file = find_iri_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );


## Checking the input from the ketchup file (e-/m²/s)
h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
μ_scatterings = setup(altitude_max, θ_lims, E_max, msis_file, iri_file);

input_type = "from_ketchup_file";
input_file = "data/Optical_Aurora_course/input_from_ketchup_3keV-10/Ie_incoming.mat"
INPUT_OPTIONS = (;input_type, input_file);
Ie_top = AURORA.Ie_top_from_ketchup(1:1:1, E, 1, μ_center, INPUT_OPTIONS.input_file)

## Checking the input from Maxwellian with LET (e-/m²/s/eV)
h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
μ_scatterings = setup(altitude_max, θ_lims, E_max, msis_file, iri_file);

input_type = "LET"
E0 = 2760;                  # characteristic energy of the maxwellian (eV)
Q = 2.02e8;                    # energy flux (eV/cm²/s)
Beams = 1:9;                # beam numbers for the precipitation, starting with field aligned down
low_energy_tail = true;     # boolean to include energy tail or not
INPUT_OPTIONS = (;input_type, E0, Q, Beams, low_energy_tail);
Ie_top = AURORA.Ie_with_LET(INPUT_OPTIONS.E0, INPUT_OPTIONS.Q, E, μ_center,
                            INPUT_OPTIONS.Beams, INPUT_OPTIONS.low_energy_tail)


## Check characteristic energy
# sum(Ie_top[1, 1, :] .* E) / sum(Ie_top[1, 1, :]) # Gaussian
sum(sum(Ie_top[1:9, 1, :]; dims=1)[1, :] .* E) / sum(Ie_top[1:9, 1, :]) # Gaussian
# Gaussian 1keV: 1280 eV
# Gaussian 3keV: 3091 eV
# IMPORTANT: Ie_with_LET flux are given in #e-/m²/s/eV
sum(sum(Ie_top[1:9, 1, :]; dims=1)[1, :] .* E .* dE) / sum(sum(Ie_top[1:9, 1, :]; dims=1)[1, :] .* dE) # Maxwellian
# Maxwellian with E0 = 1080 eV and Emax = 10 keV --> 1280 eV
# Maxwellian with E0 = 2760 eV and Emax = 20 keV --> 3090 eV

## Check total energy flux
sum(sum(Ie_top[1:9, 1, :]; dims=1)[1, :] .* E)
# Gaussian 1keV: 2.5e8
# Gaussian 3keV: 5.75e8
# IMPORTANT: Ie_with_LET flux are given in #e-/m²/s/eV
sum(sum(Ie_top[1:9, 1, :]; dims=1)[1, :] .* dE .* E)
# Maxwellian with E0 = 1080 eV, Emax = 10 keV and Q = 0.85e8 --> 2.5e8 eV
# Maxwellian with E0 = 2760 eV, Emax = 20 keV and Q = 2.02e8 --> 5.74e8 eV

##
using GLMakie
BW = AURORA.beam_weight(180:-10:0)

set_theme!(Theme(fontsize = 20, linewidth=3))
f = Figure(size = (1000, 800))
ax = Axis(f[1, 1]; xscale = log10, yscale = log10, xminorticks = IntervalsBetween(9),
          xminorticksvisible = true, xminorgridvisible = true,
          yminorgridvisible = true, ylabel = "Ie (#e-/m²/s/eV/ster)", xlabel = "E (eV)",
          title = "Maxwellian 3 keV")
lines!(E, Ie_top[1, 1, :] ./ BW[1]) # Maxwellian
# lines!(E, Ie_top[1, 1, :] ./ dE ./ BW[1]) # Gaussian
# display(f)

ax_2 = Axis(f[2, 1]; xscale = log10, yscale = log10, xminorticks = IntervalsBetween(9),
            xminorticksvisible = true, xminorgridvisible = true,
            yminorgridvisible = true, ylabel = "IeE (eV/m²/s/eV/ster)", xlabel = "E (eV)")
lines!(E, Ie_top[1, 1, :] .* E ./ BW[1]) # Maxwellian
# lines!(E, Ie_top[1, 1, :] .* E ./ dE ./ BW[1]) # Gaussian
display(GLMakie.Screen(), f)

set_theme!()


##
f, ax, l = lines(E, Ie_top[1, 1, :] ./ dE ./ BW[1])
f, ax, l = lines(E, Ie_top[1, 1, :] .* dE ./ BW[1])
lines!(E, Ie_top[1, 1, :] ./ BW[1])
ax.xscale = log10
ax.yscale = log10
f
lines!(E, Ie_top[1, 1, :] ./ BW[1] ./ π * 2)
lines(E, Ie_top[1, 1, :] ./ BW[1])
lines!(E, Ie_top[2, 1, :] ./ BW[2])
lines!(E, Ie_top[3, 1, :] ./ BW[3])
lines!(E, Ie_top[4, 1, :] ./ BW[4])
lines!(E, Ie_top[5, 1, :] ./ BW[5])
lines!(E, Ie_top[6, 1, :] ./ BW[6])
lines!(E, Ie_top[16, 1, :] ./ BW[16])
lines!(E, Ie_top[17, 1, :] ./ BW[17])
lines!(E, Ie_top[18, 1, :] ./ BW[18])
display(GLMakie.Screen(), f)
