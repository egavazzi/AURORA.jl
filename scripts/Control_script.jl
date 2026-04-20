using AURORA

## Setting parameters
altitude_lims = [100, 600];     # (km) altitude limits of the ionosphere
θ_lims = 180:-10:0              # (°) angle-limits for the electron beams
E_max = 1000;                   # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

msis_file = find_msis_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );
iri_file = find_iri_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );

## Build the model
model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)

## Define where to save the results
root_savedir = ""   # name of the root folder
name_savedir = ""   # name of the experiment folder
savedir = make_savedir(root_savedir, name_savedir)

## Define input flux
# Flickering with flat spectrum and sinusoidal modulation at 5 Hz,
# with a source set at 3000 km altitude
flux = InputFlux(FlatSpectrum(1e-2; E_min=100), SinusoidalFlickering(5.0);
                 beams=1, z_source=3000.0)

## Create and run the simulation
mode = TimeDependent(duration = 0.5,            # (s) total simulation time
                     dt = 0.001,                # (s) time step for saving data
                     CFL_number = 128,
                     # n_loop = 10,             # (optional) define manually the number of loops to run
                     # max_memory_gb = 8.0,     # (optional) or determine n_loop based on limit memory usage
                     )

sim = AuroraSimulation(model, flux, savedir; mode)
run!(sim)

## Run the analysis
make_Ie_top_file(sim)
make_volume_excitation_file(sim)
make_current_file(sim)
make_column_excitation_file(sim)
