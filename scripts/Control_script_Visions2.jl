## This is the control script from where simulations are run
using AURORA

## Setting parameters
altitude_max = 600;         # (km) top altitude of the ionosphere
θ_lims = 180:-10:0;         # (°) angle-limits for the electron beams
E_max = 3000;               # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;     # (°) angle between the B-field line and the zenith

t_sampling = 0:0.001:0.1;   # (s) time-array over which data will be saved
n_loop = 10;                # number of loops to run

Nthreads = 10;   # number of threads to be used for calculations of the energy_degradation
                # 6 threads seems to be optimal on my machine (12th Gen Intel© Core™ i7-12800HX × 16 processor)

msis_file = find_nrlmsis_file(
    year=2018, month=12, day=7, hour=11, minute=15, lat=76, lon=5, height=85:1:700
    )
iri_file = pkgdir(AURORA, "internal_data/data_electron/iri20181207.txt")   # path to the iri file

root_savedir = "Visions2"   # name for the root folder where data will be saved (data/root_savedir/)
name_savedir = "Alfven_536s_no_propagation_12eV"   # name for the actual data folder of the current experiment (data/root_savedir/name_savedir/...)
                    # if left empty, the folder will be named using the current date and time (ex: data/root_savedir/20221029-1058/...)



## Choose the input type ("from_old_matlab_file", "from_file" or "flickering")
input_type = "from_file"
input_file = "/mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_536s_correct_msis_and_scattering/Ie_incoming_536s_no_up.mat"
INPUT_OPTIONS = (;input_type, input_file);



## Run the simulation
calculate_e_transport(altitude_max, θ_lims, E_max, B_angle_to_zenith, t_sampling, n_loop,
    msis_file, iri_file, root_savedir, name_savedir, INPUT_OPTIONS, Nthreads)
