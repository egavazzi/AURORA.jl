##
using AURORA
using BenchmarkTools
using MAT

altitude_max = 600;         # (km) top altitude of the ionosphere
θ_lims = 180:-10:0;         # (°) angle-limits for the electron beams
E_max = 3000;               # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;     # (°) angle between the B-field line and the zenith

msis_file = find_nrlmsis_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );
iri_file = find_iri_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );


root_savedir = ""   # name for the root folder where data will be saved (data/root_savedir/)
name_savedir = ""   # name for the actual data folder of the current experiment (data/root_savedir/name_savedir/...)



input_type = "constant_onset"
IeE_tot = 1e-2;             # (W/m²) total energy flux of the FAB
z₀ = altitude_max;          # (km) altitude of the source
E_min = E_max - 100;         # (eV) bottom energy of the FAB
Beams = 1:2;                # beam numbers for the precipitation, starting with field aligned down
t0 = 0;                     # (s) time of start for smooth transition
t1 = 0;                     # (s) time of end for smooth transition
INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, Beams, t0, t1);


calculate_e_transport_steady_state(altitude_max, θ_lims, E_max, B_angle_to_zenith,
    msis_file, iri_file, root_savedir, name_savedir, INPUT_OPTIONS)
