# Copy of the script that was used to generate the reference results

# ======================================================================================== #
#                                  STEADY STATE                                            #
# ======================================================================================== #

using AURORA
altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
θ_lims = 180:-30:0;             # (°) angle-limits for the electron beams
E_max = 500;                   # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

msis_file = "test/regression/reference_results/msis_20051008-2200_70N-19E.txt"
iri_file = "test/regression/reference_results/iri_20051008-2200_70N-19E.txt"

## Build the model
model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)

## Define where to save the results
savedir = mktempdir()

## Define input flux
flux = InputFlux(FlatSpectrum(1e-2; E_min=E_max - 100); beams=1:2)

## Run the simulation
sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode())
run!(sim)

## Analyze the results
make_volume_excitation_file(sim)

## Overwrite the reference results
source_file = joinpath(savedir, "Qzt_all_L.mat")
dest_file = "test/regression/reference_results/SS/Qzt_all_L.mat"
cp(source_file, dest_file; force = true)







# ======================================================================================== #
##                                 TIME DEPENDENT                                          #
# ======================================================================================== #

using AURORA
altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
θ_lims = 180:-30:0;             # (°) angle-limits for the electron beams
E_max = 500;                   # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

msis_file = "test/regression/reference_results/msis_20051008-2200_70N-19E.txt"
iri_file = "test/regression/reference_results/iri_20051008-2200_70N-19E.txt"

## Build the model
model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)

## Define where to save the results
savedir = mktempdir()

## Define input flux
flux = InputFlux(FlatSpectrum(1e-2; E_min=100.0), SinusoidalFlickering(5.0);
                 beams=1, z_source=1000.0)

## Run the simulation
sim = AuroraSimulation(model, flux, savedir;
                       mode=TimeDependentMode(duration = 0.2, dt = 0.01,
                                              CFL_number = 128, n_loop = 2))
run!(sim)

## Analyze the results
make_volume_excitation_file(sim)

## Overwrite the reference results
source_file = joinpath(savedir, "Qzt_all_L.mat")
dest_file = "test/regression/reference_results/TD/Qzt_all_L.mat"
cp(source_file, dest_file; force = true)
