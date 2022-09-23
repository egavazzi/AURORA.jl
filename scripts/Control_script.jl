## This is the control script from where simulations are run
using Aurora

## Setting parameters
altitude_max = 400;         # (km) top altitude of the ionosphere
θ_lims = 180:-10:0;         # (°) angle-limits for the electron beams
E_max = 100;               # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;     # (°) angle between the B-field line and the zenith 

t = 0:0.001:0.05;           # (s) time-array over which data will be saved
n_loop = 1;                 # number of loops to run

input_file = "/mnt/data/etienne/AURORA/MI_coupling/TIME2PLAY/conversion_1.27e7-1/Ie_incoming.mat"

## Run the simulation
calculate_e_transport(altitude_max, θ_lims, E_max, B_angle_to_zenith, t, n_loop, input_file)