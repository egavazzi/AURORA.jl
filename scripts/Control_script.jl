using AURORA

## Setting parameters
altitude_lims = [100, 600];     # (km) altitude limits of the ionosphere
θ_lims = 180:-10:0              # (°) angle-limits for the electron beams
E_max = 3000;                   # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

t_sampling = 0:0.001:0.1;       # (s) time-array over which data will be saved
n_loop = 10;                    # number of loops to run

CFL_number = 128;

msis_file = find_nrlmsis_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );
iri_file = find_iri_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );

## Define where to save the results
root_savedir = ""   # name of the root folder
name_savedir = ""   # name of the experiment folder
savedir = make_savedir(root_savedir, name_savedir)

## Define input parameters
# input_type = "from_file"
# input_file = "/mnt/data/etienne/Julia/AURORA/data/MI_coupling/1.27e7-1/Ie_precipitating.mat"
# INPUT_OPTIONS = (;input_type, input_file);

input_type = "flickering";
IeE_tot = 1e-2;             # (W/m²) total energy flux of the FAB
z₀ = 3000;                  # (km) altitude of the source
E_min = 100;                # (eV) bottom energy of the FAB
f = 5;                      # (Hz) frequence of the modulation
Beams = 1;                  # beam numbers for the precipitation, starting with field aligned down
modulation = "sinus";       # type of the modulation ("square" or "sinus")
INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, f, Beams, modulation);

# input_type = "constant_onset"
# IeE_tot = 1e-2;             # (W/m²) total energy flux of the FAB
# z₀ = altitude_lims[2];      # (km) altitude of the source
# E_min = E_max - 100;        # (eV) bottom energy of the FAB
# Beams = 1:2;                # beam numbers for the precipitation, starting with field aligned down
# t0 = 0;                     # (s) time of start for smooth transition
# t1 = 0;                     # (s) time of end for smooth transition
# INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, Beams, t0, t1);

## Run the simulation
calculate_e_transport(altitude_lims, θ_lims, E_max, B_angle_to_zenith, t_sampling, n_loop,
    msis_file, iri_file, savedir, INPUT_OPTIONS, CFL_number)

## Run the analysis
make_Ie_top_file(savedir)
make_volume_excitation_file(savedir)
make_current_file(savedir)
make_column_excitation_file(savedir)
