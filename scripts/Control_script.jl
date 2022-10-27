## This is the control script from where simulations are run
using Aurora

## Setting parameters
altitude_max = 400;         # (km) top altitude of the ionosphere
θ_lims = 180:-10:0;         # (°) angle-limits for the electron beams
E_max = 1000;               # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;     # (°) angle between the B-field line and the zenith 

t = 0:0.01:1;           # (s) time-array over which data will be saved
n_loop = 1;                 # number of loops to run

root_savedir = "tests2/"

## Choose the input type ("from_file" or "flickering")
# input_type = "from_file";
# input_file = "/mnt/data/etienne/AURORA/MI_coupling/TIME2PLAY/conversion_1.27e7-1/Ie_incoming.mat";
# INPUT_OPTIONS = (;input_type, input_file);

input_type = "flickering";
IeE_tot = 1e-2;             # (W/m²) total energy flux of the FAB
z₀ = 400;                   # (km) altitude of the source
E_min = 100;                # (eV) bottom energy of the FAB
f = 10;                     # (Hz) frequence of the modulation
Beams = 1:1;                # beam numbers for the precipitation, starting with field aligned down
modulation = "sinus";       # type of the modulation ("square" or "sinus")
INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, f, Beams, modulation);

## Run the simulation
# calculate_e_transport(altitude_max, θ_lims, E_max, B_angle_to_zenith, t, n_loop, root_savedir, input_file)

calculate_e_transport(altitude_max, θ_lims, E_max, B_angle_to_zenith, t, n_loop, root_savedir, INPUT_OPTIONS)