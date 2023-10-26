## This is the control script from where simulations are run
using AURORA

## Setting parameters
altitude_max = 500;         # (km) top altitude of the ionosphere
θ_lims = 180:-10:0;         # (°) angle-limits for the electron beams
E_max = 3000;               # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;     # (°) angle between the B-field line and the zenith

t_sampling = 0:0.001:0.01;   # (s) time-array over which data will be saved
n_loop = 1;                 # number of loops to run

Nthreads = 6;   # number of threads to be used for calculations of the energy_degradation
                # 6 threads seems to be optimal on my machine (12th Gen Intel© Core™ i7-12800HX × 16 processor)

msis_file = pkgdir(AURORA, "internal_data/data_neutrals/msis20181207.txt") # path to the msis file
iri_file = pkgdir(AURORA, "internal_data/data_electron/iri20181207.txt")   # path to the iri file

root_savedir = ""   # name for the root folder where data will be saved (data/root_savedir/)
name_savedir = ""   # name for the actual data folder of the current experiment (data/root_savedir/name_savedir/...)
                    # if left empty, the folder will be named using the current date and time (ex: data/root_savedir/20221029-1058/...)

## Choose the input type ("from_old_matlab_file", "from_file" or "flickering")
# input_type = "from_old_matlab_file";
# input_file = "/mnt/data/etienne/AURORA/MI_coupling/TIME2PLAY/conversion_1.27e7-1/Ie_incoming.mat";
# INPUT_OPTIONS = (;input_type, input_file);

# input_type = "from_file"
# input_file = "/mnt/data/etienne/Julia/AURORA/data/MI_coupling/1.27e7-1/Ie_precipitating.mat"
# INPUT_OPTIONS = (;input_type, input_file);

# input_type = "flickering";
# IeE_tot = 1e-2;             # (W/m²) total energy flux of the FAB
# z₀ = 400;                   # (km) altitude of the source
# E_min = 100;                # (eV) bottom energy of the FAB
# f = 10;                     # (Hz) frequence of the modulation
# Beams = 1:2;                # beam numbers for the precipitation, starting with field aligned down
# modulation = "sinus";       # type of the modulation ("square" or "sinus")
# INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, f, Beams, modulation);

input_type = "constant_onset"
IeE_tot = 1e-2;             # (W/m²) total energy flux of the FAB
z₀ = altitude_max;          # (km) altitude of the source
E_min = E_max - 100;         # (eV) bottom energy of the FAB
Beams = 1:2;                # beam numbers for the precipitation, starting with field aligned down
t0 = 0;                     # (s) time of start for smooth transition
t1 = 0;                     # (s) time of end for smooth transition
INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, Beams, t0, t1);

## Run the simulation
calculate_e_transport(altitude_max, θ_lims, E_max, B_angle_to_zenith, t_sampling, n_loop,
    msis_file, iri_file, root_savedir, name_savedir, INPUT_OPTIONS, Nthreads)
