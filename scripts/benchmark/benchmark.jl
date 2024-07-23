using AURORA
using BenchmarkTools
using MAT

altitude_max = 600;         # (km) top altitude of the ionosphere
θ_lims = 180:-10:0;         # (°) angle-limits for the electron beams
E_max = 3000;               # (eV) upper limit to the energy grid
B_angle_to_zenith = 13;     # (°) angle between the B-field line and the zenith

t_sampling = 0:0.001:0.01;   # (s) time-array over which data will be saved
n_loop = 1;                 # number of loops to run

msis_file = find_nrlmsis_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );
iri_file = find_iri_file(
    year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
    );

input_type = "constant_onset"
IeE_tot = 1e-2;             # (W/m²) total energy flux of the FAB
z₀ = altitude_max;          # (km) altitude of the source
E_min = E_max - 100;         # (eV) bottom energy of the FAB
Beams = 1:2;                # beam numbers for the precipitation, starting with field aligned down
t0 = 0;                     # (s) time of start for smooth transition
t1 = 0;                     # (s) time of end for smooth transition
INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, Beams, t0, t1);


## Setup and initialization
h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
    μ_scatterings = setup(altitude_max, θ_lims, E_max, msis_file, iri_file);

I0 = zeros(length(h_atm) * length(μ_center), length(E));    # starting e- flux profile

t, CFL_factor = CFL_criteria(t_sampling, h_atm, v_of_E(E_max), 64)

Ie_top = Ie_top_constant(t, E, dE, n_loop, μ_center, h_atm,
                                μ_scatterings.BeamWeight, INPUT_OPTIONS.IeE_tot,
                                INPUT_OPTIONS.z₀, INPUT_OPTIONS.E_min, INPUT_OPTIONS.Beams,
                                INPUT_OPTIONS.t0, INPUT_OPTIONS.t1);

# Calculate the phase functions and put them in a Tuple
phaseN2e, phaseN2i = phase_fcn_N2(μ_scatterings.theta1, E);
phaseO2e, phaseO2i = phase_fcn_O2(μ_scatterings.theta1, E);
phaseOe, phaseOi = phase_fcn_O(μ_scatterings.theta1, E);
phase_fcn_neutrals = ((phaseN2e, phaseN2i), (phaseO2e, phaseO2i), (phaseOe, phaseOi));
cascading_neutrals = (cascading_N2, cascading_O2, cascading_O) # tuple of functions

#
i = 1

D = make_D(E, dE, θ_lims);

# Extract the top flux for the current loop
Ie_top_local = Ie_top[:, (1 + (i - 1) * (length(t) - 1)) : (length(t) + (i - 1) * (length(t) - 1)), :];

Ionization_matrix = [zeros(length(h_atm) * length(μ_center), length(t)) for _ in 1:15]
Ionizing_matrix = [zeros(length(h_atm) * length(μ_center), length(t)) for _ in 1:15]
secondary_vector = [zeros(length(E)) for _ in 1:15]
primary_vector = [zeros(length(E)) for _ in 1:15]

##
iE = length(E)
    @time A = make_A(n_neutrals, σ_neutrals, ne, Te, E, dE, iE);
    @time B, B2B_inelastic_neutrals = make_B(n_neutrals, σ_neutrals, E_levels_neutrals,
                phase_fcn_neutrals, dE, iE, μ_scatterings.Pmu2mup,
                μ_scatterings.BeamWeight_relative, μ_scatterings.theta1);

    # Compute the flux of e-
    Q  = zeros(length(h_atm) * length(μ_center), length(t), length(E));
    Ie = rand(length(h_atm) * length(μ_center), length(t), length(E));
    @btime Ie[:, :, iE] = Crank_Nicolson_Optimized(t, h_atm ./ cosd(B_angle_to_zenith), μ_center, v_of_E(E[iE]),
                                                A, B, D[iE, :], Q[:, :, iE], Ie_top_local[:, :, iE], I0[:, iE]);
        # 0.028401 seconds (3.76 k allocations: 46.271 MiB, 15.74% gc time)
        # 3000eV -> 337 energy grid points -> 337*50MiB ≈ 17GiB

    # Update the cascading of e-
    # Q  = zeros(length(h_atm) * length(μ_center), length(t), length(E));
    # @time update_Q!(Q, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals, B2B_inelastic_neutrals,
    #             cascading_neutrals, E, dE, iE, μ_scatterings.BeamWeight, μ_center)
        # 0.897071 seconds (144.07 k allocations: 228.009 MiB, 1.48% gc time)
        # ≈1.20s for highest iE, ≈220MiB for most of them (iE > 60)
        # 3000eV -> 337 energy grid points -> 300*220MiB ≈ 66GiB

    # 1st run
    # Benchmarking for some iE iterations
    Ie = rand(length(h_atm) * length(μ_center), length(t), length(E));
    Q  = zeros(length(h_atm) * length(μ_center), length(t), length(E));
    @benchmark update_Q!(Q, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals, B2B_inelastic_neutrals,
                    cascading_neutrals, E, dE, iE, μ_scatterings.BeamWeight, μ_center, 20) seconds=10
    @btime update_Q!(Q, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals, B2B_inelastic_neutrals,
                    cascading_neutrals, E, dE, iE, μ_scatterings.BeamWeight, μ_center, 20)

    using AURORA
    @benchmark update_Q_turbo!(Q, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals,
                    B2B_inelastic_neutrals, cascading_neutrals, E, dE, iE,
                    μ_scatterings.BeamWeight, μ_center,
                    Ionization_matrix, Ionizing_matrix, secondary_vector, primary_vector) seconds=20
    @time update_Q_turbo!(Q, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals,
                    B2B_inelastic_neutrals, cascading_neutrals, E, dE, iE,
                    μ_scatterings.BeamWeight, μ_center,
                    Ionization_matrix, Ionizing_matrix, secondary_vector, primary_vector)
    # b = @benchmarkable update_Q!(Q, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals, B2B_inelastic_neutrals,
    #                 cascading_neutrals, E, dE, iE, μ_scatterings.BeamWeight, μ_center)
    # run(b, seconds=30, samples=1000)
    # run(b)

    Q_save = copy(Q);
    Q ≈ Q_save
    Q[:, :, 335] .≈ Q_save[:, :, 335]


Q[:, :, 335]
Q_save[:, :, 335]

using Profile
Profile.init(n = 10^8, delay = 0.01)
    # p = Progress(length(E));
    @time for iE in length(E):-1:1
        update_Q!(Q, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals, B2B_inelastic_neutrals,
        cascading_neutrals, E, dE, iE, μ_scatterings.BeamWeight, μ_center, 20)
        # next!(p)
    end

    # p = Progress(length(E));
    @time for iE in length(E):-1:1
        update_Q_turbo!(Q, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals,
        B2B_inelastic_neutrals, cascading_neutrals, E, dE, iE,
        μ_scatterings.BeamWeight, μ_center,
        Ionization_matrix, Ionizing_matrix, secondary_vector, primary_vector)
        # next!(p)
    end

# 3keV
# size(Q) : (5742, 41, 337)
# @batch            --> 246.042693 seconds (9.41 M allocations: 22.522 GiB, 0.53% gc time)
# @turbo ionization --> 208.190795 seconds (10.04 M allocations: 22.507 GiB, 0.94% gc time)
# @turbo everywhere --> 32.799941 seconds (1.95 M allocations: 21.771 GiB, 1.30% gc time)
# size(Q) : (5742, 41, 337), no inelastic collisions
# @batch -->  31.495556 seconds (424.30 k allocations: 6.329 GiB, 0.63% gc time)
# @turbo -->  14.145960 seconds (1.11 M allocations: 6.316 GiB, 0.39% gc time)

# 5keV
# size(Q) : (5742, 81, 508)
# @batch            --> 545.570189 seconds (14.80 M allocations: 48.249 GiB, 0.36% gc time, 0.03% compilation time)
# @turbo ionization --> 385.343121 seconds (15.75 M allocations: 48.216 GiB, 0.38% gc time, 0.04% compilation time)
# @turbo everywhere -->  99.129105 seconds (2.95 M allocations: 47.063 GiB, 0.56% gc time)
# size(Q) : (5742, 81, 508), no inelastic collisions
# @batch --> 175.293155 seconds (644.11 k allocations: 18.504 GiB, 0.33% gc time)
# @turbo -->  54.258259 seconds (1.67 M allocations: 18.474 GiB, 0.80% gc time)

# 7keV
# size(Q) : (5742, 81, 680)
# @batch            --> 784.074743 seconds (19.91 M allocations: 65.325 GiB, 0.24% gc time)
# @turbo ionization --> 516.733594 seconds (21.15 M allocations: 65.266 GiB, 0.35% gc time)
# @turbo everywhere --> 158.965800 seconds (4.06 M allocations: 63.706 GiB, 0.55% gc time)
# size(Q) : (5742, 81, 680), no inelastic
# @batch --> 350.114279 seconds (952.20 k allocations: 25.082 GiB, 0.17% gc time)
# @turbo --> 106.248051 seconds (2.33 M allocations: 25.027 GiB, 0.55% gc time)



##
# Normal : 220.910 ms (47158 allocations: 90.85 MiB)
# without mul! : 168.696 ms (46358 allocations: 90.82 MiB)
# without Solve : 80.101 ms (13558 allocations: 88.74 MiB)






## Benchmark ⬇⬇⬇

# altitude_max = 600;
# θ_lims = 180:-10:0;
# E_max = 7000;
# t = 0:0.0001:0.01;

# old solve and old sparsecat : 87.272 ms (58764 allocations: 345.25 MiB)
# old solve and new sparse :    43.557 ms (10760 allocations: 58.26 MiB)
# new ldiv! + klu :     40.268 ms (8015 allocations: 54.43 MiB)
# new ldiv! + klu! :    35.911 ms (7990 allocations: 45.84 MiB)
