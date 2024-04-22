##
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PYTHONCALL_EXE"] = "/usr/local/bin/python3"  # optional
#ENV["JULIA_PYTHONCALL_EXE"] = "@PyCall"  # optional

using AURORA
# using BenchmarkTools
using MAT

msis_file = find_nrlmsis_file(
    year=2006, month=12, day=12, hour=19, minute=30, lat=70, lon=19, height=85:1:700
    );
iri_file = find_iri_file(
    year=2006, month=12, day=12, hour=19, minute=30, lat=70, lon=19, height=85:1:700
    );


root_savedir = "ionprod"   
# name for the root folder where data will be saved (data/root_savedir/)
name_savedir = "1e5keV"   
# name for the actual data folder of the current experiment (data/root_savedir/name_savedir/...)



altitude_max = 600;         # (km) top altitude of the ionosphere
θ_lims = 180:-10:0;         # (°) angle-limits for the electron beams]
#E_max = 100000;               # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;     # (°) angle between the B-field line and the zenith

input_type = "constant_onset"
IeE_tot = 1e-2;             # (W/m²) total energy flux of the FAB
z₀ = altitude_max;          # (km) altitude of the source
#E_min = E_max - 100;         # (eV) bottom energy of the FAB
Beams = 1:length(θ_lims)-1;   # beam numbers for the precipitation, starting with field aligned down
t0 = 0;                     # (s) time of start for smooth transition
t1 = 0;                     # (s) time of end for smooth transition
INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, Beams, t0, t1);



e_grid = 10 .^(range(1,stop=5,length=400))

for i in 1:length(e_grid)-1
    E_max = e_grid[i+1]
    E_min = e_grid[i]
    name_savedir = string((E_min + E_max)/2) * "eV"
    calculate_e_transport_steady_state(altitude_max, θ_lims, E_max, B_angle_to_zenith,
        msis_file, iri_file, root_savedir, name_savedir, INPUT_OPTIONS)
end