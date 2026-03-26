using AURORA

## Setting parameters
altitude_lims = [100, 600];     # (km) altitude limits of the ionosphere
θ_lims = 180:-10:0;            # (°) angle-limits for the electron beams
E_max = 3000;                   # (eV) upper limit to the energy grid
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
flux = InputFlux(FlatSpectrum(1e-2; E_min=E_max - 100); beams=1:2)

## Create and run the simulation
sim = AuroraSimulation(model, flux, savedir)
run!(sim)

## Analyze the results
make_Ie_top_file(savedir)
make_volume_excitation_file(savedir)
make_current_file(savedir)
# make_column_excitation_file(savedir) -- does not make sense for steady-state
