## This is the control script from where simulations are run
using AURORA

## Setting parameters
altitude_max = 500;         # (km) top altitude of the ionosphere
θ_lims = 180:-10:0;         # (°) angle-limits for the electron beams
E_max = 50;               # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;     # (°) angle between the B-field line and the zenith

t_sampling = 0:0.01:0.01;  # (s) time-array over which data will be saved
n_loop = 1;                 # number of loops to run

Nthreads = 20;      # number of threads to be used for calculations of the energy_degradation
                    # 6 threads seems to be optimal on my machine (12th Gen Intel© Core™ i7-12800HX × 16 processor)
CFL_number = 64;


msis_file = find_nrlmsis_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );
iri_file = find_iri_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );

# root_savedir = "CFL_simulations"   # name for the root folder where data will be saved (data/root_savedir/)
# name_savedir = "$(Int(E_max / 1000))keV_CFL-$CFL_number"   # name for the actual data folder of the current experiment (data/root_savedir/name_savedir/...)
                    # if left empty, the folder will be named using the current date and time (ex: data/root_savedir/20221029-1058/...)

root_savedir = "Benchmark"   # name for the root folder where data will be saved (data/root_savedir/)
name_savedir = "500km_$(E_max)eV_0-01(2)s"   # name for the actual data folder of the current experiment (data/root_savedir/name_savedir/...)


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
E_min = E_max - 100;        # (eV) bottom energy of the FAB
Beams = 1:2;                # beam numbers for the precipitation, starting with field aligned down
t0 = 0;                     # (s) time of start for smooth transition
t1 = 0;                     # (s) time of end for smooth transition
INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, Beams, t0, t1);



## Run the simulation
Profile.init(n = 10^7, delay = 0.01)
calculate_e_transport(altitude_max, θ_lims, E_max, B_angle_to_zenith, t_sampling, n_loop,
    msis_file, iri_file, root_savedir, name_savedir, INPUT_OPTIONS, Nthreads, CFL_number)