@testitem "Maxwellian with LET" begin
    ## Setting parameters
    altitude_lims = [100, 600];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-10:0              # (°) angle-limits for the electron beams
    E_max = 100000;                 # (eV) upper limit to the energy grid
    msis_file = find_nrlmsis_file();
    iri_file = find_iri_file();
    h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
    μ_scatterings = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file);

    # Check that IeE_top is respected for a Maxwellian without LET (with a tolerance of 0.1%)
    Ie_top = Ie_with_LET(100, 1e5, E, dE, μ_center, μ_scatterings.BeamWeight, 1, low_energy_tail = false)
    IeE_top = sum(Ie_top .* reshape(E .+ dE / 2, (1, 1, :)))
    @test isapprox(IeE_top, 1e5, rtol = 0.001)

    # Check that varying the beams does not change IeE_top
    Ie_top_LET_1 = Ie_with_LET(100, 1e5, E, dE, μ_center, μ_scatterings.BeamWeight, 1, low_energy_tail = true)
    Ie_top_LET_2 = Ie_with_LET(100, 1e5, E, dE, μ_center, μ_scatterings.BeamWeight, 1:2, low_energy_tail = true)
    Ie_top_LET_3 = Ie_with_LET(100, 1e5, E, dE, μ_center, μ_scatterings.BeamWeight, [2, 5], low_energy_tail = true)
    IeE_top_1 = sum(Ie_top_LET_1 .* reshape(E .+ dE / 2, (1, 1, :)))
    IeE_top_2 = sum(Ie_top_LET_2 .* reshape(E .+ dE / 2, (1, 1, :)))
    IeE_top_3 = sum(Ie_top_LET_3 .* reshape(E .+ dE / 2, (1, 1, :)))

    @test IeE_top_1 ≈ IeE_top_2 ≈ IeE_top_3
end


@testitem "(SS) Does LET run?" begin
    # Setting parameters
    altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-45:0;            # (°) angle-limits for the electron beams
    E_max = 100;                   # (eV) upper limit to the energy grid
    B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith
    msis_file = find_nrlmsis_file();
    iri_file = find_iri_file();
    savedir = make_savedir("", "")
    # Define input parameters
    input_type = "LET"
    E0 = 50                 # characteristic energy (eV)
    Q = 6.24151e15          # total energy flux (eV/m2/s) (1 mW = 24151e15 eV)
    Beams = 1               # indices of the beams with precipitation
    low_energy_tail = true  # with or without a low energy tail
    INPUT_OPTIONS = (;input_type, E0, Q, Beams, low_energy_tail);
    # Run simulation
    calculate_e_transport_steady_state(altitude_lims, θ_lims, E_max, B_angle_to_zenith,
                                    msis_file, iri_file, savedir, INPUT_OPTIONS)
    @test true
end

@testitem "(SS) Does constant onset run?" begin
    # Setting parameters
    altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-45:0;            # (°) angle-limits for the electron beams
    E_max = 100;                 # (eV) upper limit to the energy grid
    B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith
    msis_file = find_nrlmsis_file();
    iri_file = find_iri_file();
    savedir = make_savedir("", "")
    # Define input parameters
    input_type = "constant_onset"
    IeE_tot = 1.0           # total energy flux (W/m²)
    z₀ = 500.0              # altitude of the source (km)
    E_min = 50.0            # minimum energy (eV)
    Beams = 1:2             # indices of the beams with precipitation
    t0 = 0.0                # onset start time (s)
    t1 = 0.0                # onset end time (s)
    INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, Beams, t0, t1);
    # Run simulation
    calculate_e_transport_steady_state(altitude_lims, θ_lims, E_max, B_angle_to_zenith,
                                    msis_file, iri_file, savedir, INPUT_OPTIONS)
    @test true
end

@testitem "(TD) Does constant onset run?" begin
    # Setting parameters
    altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-45:0;            # (°) angle-limits for the electron beams
    E_max = 100;                 # (eV) upper limit to the energy grid
    B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith
    t_sampling = 0:0.01:0.1;       # (s) time-array over which data will be saved
    n_loop = 2;                    # number of loops to run
    CFL_number = 128;
    msis_file = find_nrlmsis_file();
    iri_file = find_iri_file();
    savedir = make_savedir("", "")
    # Define input parameters
    input_type = "constant_onset"
    IeE_tot = 1.0           # total energy flux (W/m²)
    z₀ = 500.0              # altitude of the source (km)
    E_min = 50.0            # minimum energy (eV)
    Beams = 1:2             # indices of the beams with precipitation
    t0 = 0.0                # onset start time (s)
    t1 = 0.05               # onset end time (s)
    INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, Beams, t0, t1);
    # Run simulation
    calculate_e_transport(altitude_lims, θ_lims, E_max, B_angle_to_zenith, t_sampling, n_loop,
                          msis_file, iri_file, savedir, INPUT_OPTIONS, CFL_number)
    @test true
end

@testitem "(TD) Does flickering run?" begin
    # Setting parameters
    altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-45:0;            # (°) angle-limits for the electron beams
    E_max = 100;                 # (eV) upper limit to the energy grid
    B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith
    t_sampling = 0:0.01:0.1;       # (s) time-array over which data will be saved
    n_loop = 2;                    # number of loops to run
    CFL_number = 128;
    msis_file = find_nrlmsis_file();
    iri_file = find_iri_file();
    savedir = make_savedir("", "")
    # Define input parameters
    input_type = "flickering";
    IeE_tot = 1e-2;             # (W/m²) total energy flux of the FAB
    z₀ = 1000;                  # (km) altitude of the source
    E_min = 50;                # (eV) bottom energy of the FAB
    f = 5;                      # (Hz) frequence of the modulation
    Beams = 1;                  # beam numbers for the precipitation, starting with field aligned down
    modulation = "sinus";       # type of the modulation ("square" or "sinus")
    INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, f, Beams, modulation);
    # Run simulation
    calculate_e_transport(altitude_lims, θ_lims, E_max, B_angle_to_zenith, t_sampling, n_loop,
                          msis_file, iri_file, savedir, INPUT_OPTIONS, CFL_number)
    @test true
end




## Code to reproduce Fig 4. of Meier et al. 1989
#=
using CairoMakie

"""
# Note
For things to work with the factors given in b, everything should be in keV

# Arguments
- E: energy (keV)
- Q: total energy flux (keV/cm²/s)
- E0: characteristic energy (keV)

# Returns
- Ie: total primary electron spectra
"""
function equations_from_paper(E, E0, Q; low_energy_tail = true)
    # Parameter for LET amplitude (units of keV)
    b = 1 / (0.8 * E0) .* (E0 < 0.5) +
        1 / (0.1 * E0 + 0.35) .* (E0 >= 0.5)

    # Maxwellian spectra (e⁻/cm²/s/keV/ster)
    Ie = Q / (2π * E0^3) .* E .* exp.(-E ./ E0)

    if low_energy_tail
        Ie_max = maximum(Ie)
        Ie = Ie .+ 0.4 * Ie_max * (E0 ./ E) .* exp.(-E ./ b)
    end

    return Ie # in (e⁻/cm²/s/keV/ster)
end

function equations_modified(E, E0, Q; low_energy_tail = true)
    # Parameter for LET amplitude (units of keV)
    b = (0.8 * E0) .* (E0 < 0.5) +
        (0.1 * E0 + 0.35) .* (E0 >= 0.5)

    # Maxwellian spectra (e⁻/cm²/s/keV/ster)
    Ie = Q / (2π * E0^3) .* E .* exp.(-E ./ E0)

    if low_energy_tail
        Ie_max = maximum(Ie)
        Ie = Ie .+ 0.4 * Ie_max * (E0 ./ E) .* exp.(-E ./ b)
    end

    return Ie # in (e⁻/cm²/s/keV/ster)
end




dE = 0.01 # energy grid step size (keV)
E = 0.01:dE:1e2 # energy grid (keV)
E0_to_run = [0.1, 0.25, 0.50, 1, 2, 4, 8] # characteristic energies (keV)
Q = 6.242e8 # total downard energy flux of 1.0 erg/cm²/s, but in keV/cm²/s

# Plot equations_from_paper
f = Figure()
ax = Axis(f[1, 1];
          xscale = log10,
          yscale = log10,
          xlabel = "Energy (eV)",
          ylabel = "Ie (e⁻/cm²/s/eV/ster)",
          title = "Equations not modified",
          )
for (i, E0) in enumerate(E0_to_run)
    Ie_top = equations_from_paper(E, E0, Q; low_energy_tail = false)
    Ie_top_LET = equations_from_paper(E, E0, Q; low_energy_tail = true)
    # We plot E in (eV) and Ie in (e⁻/cm²/s/eV/ster)
    lines!(E * 1e3, Ie_top / 1e3 * pi, color = Cycled(i), linestyle = :dash)
    lines!(E * 1e3, Ie_top_LET / 1e3 * pi, color = Cycled(i))
    # Sanity check:
    # Calculate the total enery flux in (keV/cm²/s) from the flux in (e⁻/cm²/s/keV/ster)
    # We need to multiply by E and by dE
    # Should be equal to Q for the simple Maxwellian without LET
    println(sum(Ie_top .* dE .* E  * pi))
end
#=
Their normalization in ster looks a bit weird to me.
In order to find Q back, I had to multiply the flux by π instead of 2π (half a sphere).
Also, the flux in fig4 do not seem to be in /ster, as I had to multiply them by π to
get similar values as the figure.
Anyway, this just shifts the plot slightly.
=#
ylims!(10, 1e8)
display(f)


# Plot equations_modified
f = Figure()
ax = Axis(f[1, 1];
          xscale = log10,
          yscale = log10,
          xlabel = "Energy (eV)",
          ylabel = "Ie (e⁻/cm²/s/eV/ster)",
          title = "Equations modified",
          )
for (i, E0) in enumerate(E0_to_run)
    Ie_top = equations_modified(E, E0, Q; low_energy_tail = false)
    Ie_top_LET = equations_modified(E, E0, Q; low_energy_tail = true)
    # We plot E in (eV) and Ie in (e⁻/cm²/s/eV/ster)
    lines!(E * 1e3, Ie_top / 1e3 * pi, color = Cycled(i), linestyle = :dash)
    lines!(E * 1e3, Ie_top_LET / 1e3 * pi, color = Cycled(i))
end
ylims!(10, 1e8)
display(f)
=#
