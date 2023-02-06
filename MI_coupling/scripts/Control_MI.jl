"""
This script controls the MI coupling.
It calls Aurora.jl, ketchup, and some MatLab functions from AURORA
"""

using MATLAB
using Aurora

current_folder = pwd()

# Path to Vlasov simulation
path_to_vlasov_initial_file = "/mnt/data/etienne/ketchup/MI_coupling/l7/1.27e7/Nz_3000_420to480s/outp/fzvzmuIB2400000.mat"
path_to_vlasov_simulation = "/mnt/data/etienne/ketchup/MI_coupling/l7/1.27e7/Nz_3000_420to480s/"
path_to_next_vlasov_simulation = "/mnt/data/etienne/ketchup/MI_coupling/l7/1.27e7_MI/0.3s-1"
# Path to Aurora simulation
path_to_aurora_simulation = "/mnt/data/etienne/Julia/Aurora/data/MI_coupling/1.27e7-2s_1/"


## 1. Convert .dat ketchup results to .mat
# cd(path_to_vlasov_simulation)
@mput path_to_vlasov_simulation
mat"cd path_to_vlasov_simulation"
mat"ketchup_b6conv" # run ketchup_b6conv.m in a Matlab engine
# cd(current_folder)

## 2. Convert fzvzmu from ketchup into Ie that can be used in Aurora.jl
index_specie = 1
HMR_VZ = 20         # How Much do you want to Reduce vz
HMR_MU = 20         # How Much do you want to Reduce μ_mag
E_max = 7000        # (eV) maximum energy for the E-grid
θ_lims = 180:-10:0  # pitch-angle limits for the streams (180° is straight downward)
first_run = 1       # if 1 it will extract only the initial state (constant incoming flux)

Ie = convert_M_to_I(path_to_vlasov_initial_file, path_to_vlasov_simulation,
                        index_specie, HMR_VZ, HMR_MU, E_max, θ_lims, first_run);

cd(path_to_vlasov_simulation)
    file = matopen("Ie_precipitating.mat", "w")
        write(file, "Ie_total", Ie)
    close(file)
cd(path_to_aurora_simulation)
    file = matopen("Ie_precipitating.mat", "w")
        write(file, "Ie_total", Ie)
    close(file)
    file = open("linked_vlasov_simulations.txt", "w")
        write(file, "initial = $path_to_vlasov_initial_file \n")
        write(file, "folder_from = $path_to_vlasov_simulation \n")
    close(file)
cd(current_folder)

## 3. Calculate the ionospheric response using Aurora.jl
altitude_max = 400;         # (km) top altitude of the ionosphere
θ_lims = 180:-10:0;         # (°) angle-limits for the electron beams
E_max = 7000;               # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;     # (°) angle between the B-field line and the zenith

t = 0:0.001:0.2;            # (s) time-array over which data will be saved
n_loop = 10;                # number of loops to run

root_savedir = "MI_coupling"
name_savedir = "1.27e7-2s_1"

input_type = "from_file"
input_file = joinpath(path_to_aurora_simulation, "Ie_precipitating.mat")
INPUT_OPTIONS = (;input_type, input_file);

calculate_e_transport(altitude_max, θ_lims, E_max, B_angle_to_zenith, t, n_loop, root_savedir,
                        name_savedir, INPUT_OPTIONS)

## 4. Extract Ie at the top of the ionosphere
mat"cd /mnt/data/etienne/Julia/Aurora/data/MI_coupling/"
mat"current_folder = pwd"
mat"make_all_Ie_top_MI"

## 5. Prepare the next Vlasov simulation
if !isdir(path_to_next_vlasov_simulation)
    println("Creating folder for the next Vlasov simulation...")
    command = Cmd(`cp -r $path_to_vlasov_simulation $path_to_next_vlasov_simulation`)
    run(command)
    cd(path_to_next_vlasov_simulation)
    run(`rm dumps/*`)
    run(`find outp/ -name "*.mat" -delete`)
    run(`find outp/ -name "*.dat" -delete`)
    println("... created!")
    cd(current_folder)
end

## 6. Convert Ie from Aurora.jl to fzvzmu that can be used in ketchup

index_specie = 1
HMR_MU = 20
HMR_E = 5
E_max = 7000
θ_lims = 180:-10:0

f1, f2 = convert_I_to_M(path_to_vlasov_simulation, path_to_aurora_simulation,
                        index_specie, HMR_MU, HMR_E, E_max, θ_lims);

write(joinpath(path_to_aurora_simulation, "f_response.bin"), f1)
write(joinpath(path_to_aurora_simulation, "f_response2.bin"), f2)
write(joinpath(path_to_next_vlasov_simulation, "f_response.bin"), f1)
write(joinpath(path_to_next_vlasov_simulation, "f_response2.bin"), f2)
































cd(path_to_aurora_simulation)
    file = open("linked_vlasov_simulations.txt", "a")
            write(file, "folder_to = $path_to_next_vlasov_simulation")
    close(file)
cd(current_folder)
