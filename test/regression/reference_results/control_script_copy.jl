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
root_savedir = "temp_results/"   # name of the root folder
name_savedir = "SS/"   # name of the experiment folder
savedir = make_savedir(root_savedir, name_savedir; behavior = "custom")

## Define input flux
flux = InputFlux(FlatSpectrum(1e-2; E_min=E_max - 100); beams=1:2)

## Run the simulation
sim = AuroraSimulation(model, flux, savedir; solver=SteadyStateSolver())
run!(sim)

## Analyze the results
make_volume_excitation_file(sim)
make_column_excitation_file(sim)
make_current_file(sim)










# ======================================================================================== #
##                                 TIME DEPENDENT                                          #
# ======================================================================================== #

using AURORA
altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
θ_lims = 180:-30:0;             # (°) angle-limits for the electron beams
E_max = 500;                   # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

t_total = 0.2;                  # (s) total simulation time
dt = 0.01;                      # (s) time step for saving data
CFL_number = 128;

msis_file = "test/regression/reference_results/msis_20051008-2200_70N-19E.txt"
iri_file = "test/regression/reference_results/iri_20051008-2200_70N-19E.txt"

## Build the model
model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)

## Define where to save the results
root_savedir = "temp_results/"   # name of the root folder
name_savedir = "TD/"   # name of the experiment folder
savedir = make_savedir(root_savedir, name_savedir; behavior = "custom")

## Define input flux
flux = InputFlux(FlatSpectrum(1e-2; E_min=100.0), SinusoidalFlickering(5.0);
                 beams=1, z_source=1000.0)

## Run the simulation
sim = AuroraSimulation(model, flux, savedir;
                       solver=TimeDependentSolver(t_total, dt; CFL_number, n_loop=2))
run!(sim)

## Analyze the results
make_volume_excitation_file(sim)
make_column_excitation_file(sim)
make_current_file(sim)
