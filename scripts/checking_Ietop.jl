using AURORA

## Setting parameters
altitude_max = 600;        # (km) top altitude of the ionosphere
θ_lims = [180, 170, 150, 120, 100, 90, 80, 60, 30, 10, 0];         # (°) angle-limits for the electron beams
E_max = 8000;               # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;     # (°) angle between the B-field line and the zenith

t_sampling = 0:0.001:0.025;  # (s) time-array over which data will be saved
n_loop = 100;                 # number of loops to run

CFL_number = 128;


msis_file = find_nrlmsis_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );
iri_file = find_iri_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );


##
h_atm, ne, Te, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
μ_scatterings = setup_new(altitude_max, θ_lims, E_max, msis_file, iri_file);

I0 = zeros(length(h_atm) * length(μ_center), length(E));    # starting e- flux profile

t, CFL_factor = CFL_criteria(t_sampling, h_atm, v_of_E(E_max), CFL_number)


Ietop_file = "/mnt/data/bjorn/Campaigns/Pulsating-aurora/IeMaxTopAtm-02.mat"
Ie_top = Ie_top_from_file_Silje(t, E, μ_center, n_loop, Ietop_file)

size(Ie_top)

##
using CairoMakie
E, _ = AURORA.make_energy_grid(8000)
f, ax, hm = heatmap(range(0, 2.5, 7501), E, Ie_top[1, :, :])
