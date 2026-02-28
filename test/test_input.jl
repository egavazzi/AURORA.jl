@testitem "Ie_top_from_file" begin
    using MAT: matwrite

    # Dummy grids and parameters
    n_μ = 2
    n_E = 4
    n_t_file = 5
    μ_center = [-0.8, 0.8] # beam 1 is down-going, beam 2 is up-going
    E = [100, 200, 300, 400]

    # Helper: write a .mat file and return its path
    function write_mat(filename, Ie; t_top=nothing)
        data = Dict{String,Any}("Ie_total" => Ie)
        if !isnothing(t_top)
            data["t_top"] = t_top
        end
        matwrite(filename, data)
        return filename
    end

    # Reference flux
    ref_flux = ones(n_μ, n_t_file, n_E)

    mktempdir() do tmpdir
        ## ── Test 1: constant flux (no t_top in file) ────────────────────────────────
        path = write_mat(joinpath(tmpdir, "test_constant.mat"),
                        ones(n_μ, 1, n_E))
        t_sim = 0.0:0.01:0.1
        Ie = Ie_top_from_file(t_sim, E, μ_center, path)
        @test size(Ie) == (n_μ, length(t_sim), n_E)
        @test all(Ie[1, :, :] .== 1.0)      # down-going beam kept
        @test all(Ie[2, :, :] .== 0.0)      # up-going beam zeroed

        ## ── Test 2: matching time grids (same dt, same length) ──────────────────────
        t_file = collect(0.0:0.01:0.04)
        path = write_mat(joinpath(tmpdir, "test_match.mat"), ref_flux; t_top=t_file)
        Ie = Ie_top_from_file(t_file, E, μ_center, path) # use same time grid
        @test size(Ie) == (n_μ, length(t_file), n_E)
        @test all(Ie[1, :, :] .== 1.0)
        @test all(Ie[2, :, :] .== 0.0)

        ## ── Test 3: file dt coarser than simulation dt ───────────────────────────────
        # file: 5 steps with dt=0.02 s  →  simulation: 9 steps with dt=0.01 s
        # (non-integer ratio tested in Test 5)
        t_file = collect(0.0:0.02:0.08)
        path = write_mat(joinpath(tmpdir, "test_coarse.mat"), ref_flux; t_top=t_file)
        t_sim = 0.0:0.01:0.08
        Ie = Ie_top_from_file(t_sim, E, μ_center, path)
        @test size(Ie) == (n_μ, length(t_sim), n_E)
        @test all(Ie[1, :, :] .== 1.0)
        @test all(Ie[2, :, :] .== 0.0)

        ## ── Test 4: file dt finer than simulation dt ────────────────────────────────
        # file: 5 steps with dt=0.005 s  →  simulation: 3 steps with dt=0.01 s
        n_t_fine = 5
        t_file = collect(0.0:0.005:0.02)
        flux_fine = ones(n_μ, n_t_fine, n_E)
        flux_fine[2, :, :] .= 0.0
        path = write_mat(joinpath(tmpdir, "test_fine.mat"), flux_fine; t_top=t_file)
        t_sim = 0.0:0.01:0.02
        # A warning should be issued because file dt < sim dt
        Ie = @test_warn r"finer" Ie_top_from_file(t_sim, E, μ_center, path)
        @test size(Ie) == (n_μ, length(t_sim), n_E)
        @test all(Ie[1, :, :] .== 1.0)
        @test all(Ie[2, :, :] .== 0.0)

        ## ── Test 5: non-integer dt ratio ─────────────────────────────────────────────
        # file dt = 0.03 s, simulation dt = 0.01 s  → ratio = 3
        # file dt = 0.03 s, simulation dt = 0.02 s  → ratio = 1.5 (non-integer)
        n_t_nonint = 4
        t_file = collect(0.0:0.03:0.09)
        flux_nonint = ones(n_μ, n_t_nonint, n_E)
        flux_nonint[2, :, :] .= 0.0
        path = write_mat(joinpath(tmpdir, "test_nonint.mat"), flux_nonint; t_top=t_file)
        t_sim = 0.0:0.02:0.08
        Ie = Ie_top_from_file(t_sim, E, μ_center, path)
        @test size(Ie) == (n_μ, length(t_sim), n_E)
        @test all(Ie[1, :, :] .== 1.0)

        ## ── Test 6: file shorter than simulation → zero padding + warning ────────────
        t_file = collect(0.0:0.01:0.02)     # 3 steps, ends at 0.02 s
        n_t_short = 3
        flux_short = ones(n_μ, n_t_short, n_E)
        flux_short[2, :, :] .= 0.0
        path = write_mat(joinpath(tmpdir, "test_short.mat"), flux_short; t_top=t_file)
        t_sim = 0.0:0.01:0.09               # 10 steps, ends at 0.09 s  (longer than file)
        Ie = @test_warn r"shorter" Ie_top_from_file(t_sim, E, μ_center, path)
        @test size(Ie) == (n_μ, length(t_sim), n_E)
        @test all(Ie[1, 1:3, :] .== 1.0) # first 3 steps should have flux = 1
        @test all(Ie[1, 4:end, :] .== 0.0) # but the rest should have flux = 0

        ## ── Test 7: file longer than simulation → truncated ─────────────────────────
        t_file = collect(0.0:0.01:0.04)     # 5 steps
        path = write_mat(joinpath(tmpdir, "test_long.mat"), ref_flux; t_top=t_file)
        t_sim = 0.0:0.01:0.02               # only 3 steps
        Ie = Ie_top_from_file(t_sim, E, μ_center, path)
        @test size(Ie) == (n_μ, length(t_sim), n_E)
        @test all(Ie[1, :, :] .== 1.0)

        ## ── Test 8: interpolation keywords ──────────────────────────────────────────
        # Ramp flux from 0 to 1 over 5 file steps; check that linear interp gives
        # intermediate values on a finer simulation grid
        t_file = collect(0.0:0.1:0.4)       # 5 steps
        flux_ramp = zeros(n_μ, 5, n_E)
        for k in 1:5
            flux_ramp[1, k, :] .= (k - 1) / 4.0   # 0, 0.25, 0.5, 0.75, 1.0
        end
        path = write_mat(joinpath(tmpdir, "test_linear.mat"), flux_ramp; t_top=t_file)
        t_sim = 0.0:0.05:0.4                # finer grid, 9 steps
        Ie_const  = Ie_top_from_file(t_sim, E, μ_center, path; interpolation=:constant)
        Ie_linear = Ie_top_from_file(t_sim, E, μ_center, path; interpolation=:linear)
        # At t=0.05 (between 0.0 and 0.1 in file):
        #   constant holds 0.0  →  Ie_const[1, 2, 1] ≈ 0.0
        #   linear gives 0.125  →  Ie_linear[1, 2, 1] ≈ 0.125
        @test Ie_const[1, 2, 1] ≈ 0.0
        @test Ie_linear[1, 2, 1] ≈ 0.125

        ## ── Test 9: energy dim in file larger than simulation E grid → trimmed ──────
        n_E_big = n_E + 3
        flux_big = ones(n_μ, n_t_file, n_E_big)
        path = write_mat(joinpath(tmpdir, "test_bigE.mat"), flux_big; t_top=collect(0.0:0.01:0.04))
        Ie = Ie_top_from_file(0.0:0.01:0.04, E, μ_center, path)
        @test size(Ie, 3) == n_E
        @test all(Ie[1, :, :] .== 1.0)
        @test all(Ie[2, :, :] .== 0.0)
    end # of mktempdir
end

@testitem "Maxwellian with LET" begin
    ## Setting parameters
    altitude_lims = [100, 600];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-10:0              # (°) angle-limits for the electron beams
    E_max = 100000;                 # (eV) upper limit to the energy grid
    msis_file = find_msis_file();
    iri_file = find_iri_file();
    h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
    μ_scatterings = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file);

    # Physical constants
    qₑ = 1.602176620898e-19  # Elementary charge (C)

    # Check that IeE_top is respected for a Maxwellian without LET (with a tolerance of 0.1%)
    IeE_tot = 1e-2  # W/m²
    Ie_top = Ie_with_LET(IeE_tot, 100, E, dE, μ_center, μ_scatterings.BeamWeight, 1, low_energy_tail = false)
    IeE_top_check = sum(Ie_top .* reshape(E .+ dE / 2, (1, 1, :))) * qₑ
    @test isapprox(IeE_top_check, IeE_tot, rtol = 0.001)

    # Check that varying the beams does not change IeE_top
    Ie_top_LET_1 = Ie_with_LET(IeE_tot, 100, E, dE, μ_center, μ_scatterings.BeamWeight, 1, low_energy_tail = true)
    Ie_top_LET_2 = Ie_with_LET(IeE_tot, 100, E, dE, μ_center, μ_scatterings.BeamWeight, 1:2, low_energy_tail = true)
    Ie_top_LET_3 = Ie_with_LET(IeE_tot, 100, E, dE, μ_center, μ_scatterings.BeamWeight, [2, 5], low_energy_tail = true)
    IeE_top_1 = sum(Ie_top_LET_1 .* reshape(E .+ dE / 2, (1, 1, :))) * qₑ
    IeE_top_2 = sum(Ie_top_LET_2 .* reshape(E .+ dE / 2, (1, 1, :))) * qₑ
    IeE_top_3 = sum(Ie_top_LET_3 .* reshape(E .+ dE / 2, (1, 1, :))) * qₑ

    @test IeE_top_1 ≈ IeE_top_2 ≈ IeE_top_3
end

@testitem "Ie_top_modulated - flat spectrum" begin
    ## Setting parameters
    altitude_lims = [100, 600];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-10:0              # (°) angle-limits for the electron beams
    E_max = 5000;                 # (eV) upper limit to the energy grid
    msis_file = find_msis_file();
    iri_file = find_iri_file();
    h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
    μ_scatterings = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file);

    # Physical constants
    qₑ = 1.602176620898e-19  # Elementary charge (C)

    # Check that IeE_top is respected for different cases
    IeE_tot = 1e-2  # W/m²
    E_min = E_max - 100
    Beams = [1, 2]
    # SS
    Ie_top = Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, μ_scatterings.BeamWeight, 1:1:1, 1, h_atm;
                              spectrum=:flat, E_min=E_min)
    IeE_top_check = sum(Ie_top .* reshape(E .+ dE / 2, (1, 1, :))) * qₑ
    @test isapprox(IeE_top_check, IeE_tot, rtol = 0.001)
    # SS, different Beams
    Ie_top = Ie_top_modulated(IeE_tot, E, dE, μ_center, [1, 3, 4], μ_scatterings.BeamWeight, 1:1:1, 1, h_atm;
                              spectrum=:flat, E_min=E_min)
    IeE_top_check = sum(Ie_top .* reshape(E .+ dE / 2, (1, 1, :))) * qₑ
    @test isapprox(IeE_top_check, IeE_tot, rtol = 0.001)
    # SS, different E_min
    Ie_top = Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, μ_scatterings.BeamWeight, 1:1:1, 1, h_atm;
                              spectrum=:flat, E_min=E_max - 1000)
    IeE_top_check = sum(Ie_top .* reshape(E .+ dE / 2, (1, 1, :))) * qₑ
    @test isapprox(IeE_top_check, IeE_tot, rtol = 0.001)
    # SS, different E_min
    Ie_top = Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, μ_scatterings.BeamWeight, 1:1:1, 1, h_atm;
                              spectrum=:flat, E_min=E_max - 10)
    IeE_top_check = sum(Ie_top .* reshape(E .+ dE / 2, (1, 1, :))) * qₑ
    @test isapprox(IeE_top_check, IeE_tot, rtol = 0.001)
    # SS, different E_min
    Ie_top = Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, μ_scatterings.BeamWeight, 1:1:1, 1, h_atm;
                              spectrum=:flat, E_min=E_max)
    IeE_top_check = sum(Ie_top .* reshape(E .+ dE / 2, (1, 1, :))) * qₑ
    @test isapprox(IeE_top_check, IeE_tot, rtol = 0.001)
    # TD
    Ie_top = Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, μ_scatterings.BeamWeight, 0:0.01:1, 1, h_atm;
                              spectrum=:flat, E_min=E_min)
    IeE_top_check = sum(Ie_top[:, 1, :] .* reshape(E .+ dE / 2, (1, :))) * qₑ
    @test isapprox(IeE_top_check, IeE_tot, rtol = 0.001)
    # TD, high source altitude
    Ie_top = Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, μ_scatterings.BeamWeight, 0:0.01:1, 1, h_atm;
                              spectrum=:flat, E_min=E_min, z_source=1000.0)
    IeE_top_check = sum(Ie_top[:, end, :] .* reshape(E .+ dE / 2, (1, :))) * qₑ
    @test isapprox(IeE_top_check, IeE_tot, rtol = 0.001)
    # TD, delayed smooth onset
    Ie_top = Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, μ_scatterings.BeamWeight, 0:0.01:1, 1, h_atm;
                              spectrum=:flat, E_min=E_min, z_source=altitude_lims[2],
                              t_start=0.4, t_end=0.8)
    IeE_top_check = sum(Ie_top[:, end, :] .* reshape(E .+ dE / 2, (1, :))) * qₑ
    @test isapprox(IeE_top_check, IeE_tot, rtol = 0.001)
end


@testitem "Ie_top_modulated - gaussian spectrum" begin
    ## Setting parameters
    altitude_lims = [100, 600];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-10:0              # (°) angle-limits for the electron beams
    E_max = 10000;                  # (eV) upper limit to the energy grid
    msis_file = find_msis_file();
    iri_file = find_iri_file();
    h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
    μ_scatterings = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file);

    # Physical constants
    qₑ = 1.602176620898e-19  # Elementary charge (C)

    IeE_tot = 1e-2  # W/m²
    E₀ = 5000.0     # center energy (eV)
    ΔE = 500.0      # energy width (eV)
    Beams = [1, 2]

    # SS - check energy flux is conserved
    Ie_top = Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, μ_scatterings.BeamWeight, 1:1:1, 1, h_atm;
                              spectrum=:gaussian, E₀=E₀, ΔE=ΔE)
    IeE_top_check = sum(Ie_top .* reshape(E .+ dE / 2, (1, 1, :))) * qₑ
    @test isapprox(IeE_top_check, IeE_tot, rtol = 0.001)

    # SS, different Beams
    Ie_top = Ie_top_modulated(IeE_tot, E, dE, μ_center, [1, 3, 4], μ_scatterings.BeamWeight, 1:1:1, 1, h_atm;
                              spectrum=:gaussian, E₀=E₀, ΔE=ΔE)
    IeE_top_check = sum(Ie_top .* reshape(E .+ dE / 2, (1, 1, :))) * qₑ
    @test isapprox(IeE_top_check, IeE_tot, rtol = 0.001)

    # SS, different E₀ and ΔE
    Ie_top = Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, μ_scatterings.BeamWeight, 1:1:1, 1, h_atm;
                              spectrum=:gaussian, E₀=3000.0, ΔE=1000.0)
    IeE_top_check = sum(Ie_top .* reshape(E .+ dE / 2, (1, 1, :))) * qₑ
    @test isapprox(IeE_top_check, IeE_tot, rtol = 0.001)

    # TD
    Ie_top = Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, μ_scatterings.BeamWeight, 0:0.01:1, 1, h_atm;
                              spectrum=:gaussian, E₀=E₀, ΔE=ΔE)
    IeE_top_check = sum(Ie_top[:, 1, :] .* reshape(E .+ dE / 2, (1, :))) * qₑ
    @test isapprox(IeE_top_check, IeE_tot, rtol = 0.001)
end


@testitem "Ie_top_modulated - modulation types" begin
    ## Setting parameters
    altitude_lims = [100, 600];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-10:0              # (°) angle-limits for the electron beams
    E_max = 5000;                   # (eV) upper limit to the energy grid
    msis_file = find_msis_file();
    iri_file = find_iri_file();
    h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
    μ_scatterings = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file);

    # Physical constants
    qₑ = 1.602176620898e-19  # Elementary charge (C)

    IeE_tot = 1e-2  # W/m²
    E_min = 4000.0
    Beams = [1, 2]
    f = 10.0  # Hz

    # Sinus modulation
    Ie_top = Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, μ_scatterings.BeamWeight, 0:0.001:0.5, 1, h_atm;
                              spectrum=:flat, E_min=E_min, modulation=:sinus, f=f, amplitude=0.8)
    # At t where modulation is at peak (sin²(πft) = 1), flux should be IeE_tot
    # This happens when πft = π/2, i.e., t = 0.05s for f=10Hz
    t_peak_idx = round(Int, 0.05 / 0.001) + 1  # index for t=0.05s
    IeE_at_peak = sum(Ie_top[:, t_peak_idx, :] .* reshape(E .+ dE / 2, (1, :))) * qₑ
    @test isapprox(IeE_at_peak, IeE_tot, rtol = 0.001)
    # At t where modulation is at minimum (sin²(πft) = 0), flux should be (1 - amplitude)*IeE_tot
    # This happens when πft = 0, i.e. t = 0.1s for f=10Hz
    t_min_idx = round(Int, 0.1 / 0.001) + 1  # index for t=0.1s
    IeE_at_min = sum(Ie_top[:, t_min_idx, :] .* reshape(E .+ dE / 2, (1, :))) * qₑ
    @test isapprox(IeE_at_min, 0.2 * IeE_tot, rtol = 0.001)

    # Square modulation
    Ie_top = Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, μ_scatterings.BeamWeight, 0:0.001:0.5, 1, h_atm;
                              spectrum=:flat, E_min=E_min, modulation=:square, f=f, amplitude=0.2)
    IeE_in_time = sum(Ie_top .* reshape(E .+ dE / 2, (1, 1, :)) * qₑ, dims = (1, 3))
    # The maximum value should be equal to IeE_tot
    @test isapprox(maximum(IeE_in_time), IeE_tot, rtol=0.001)
    # The minimum value should be equal to 0.8 * IeE_tot
    @test isapprox(minimum(IeE_in_time), 0.8 * IeE_tot, rtol=0.001)

    # Partial amplitude modulation (amplitude=0.5)
    Ie_top_partial = Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, μ_scatterings.BeamWeight, 0:0.001:0.5, 1, h_atm;
                                      spectrum=:flat, E_min=E_min, modulation=:sinus, f=f, amplitude=0.5)
    # At minimum of sinus (where sin²(πft)=0), flux should be 0.5*IeE_tot
    t_min_idx = 1 # at t=0, sin²(0)=0 (modulation minimum)
    IeE_at_min = sum(Ie_top_partial[:, t_min_idx, :] .* reshape(E .+ dE / 2, (1, :))) * qₑ
    @test isapprox(IeE_at_min, 0.5 * IeE_tot, rtol=0.001)
end


@testitem "(SS) Does LET run?" begin
    # Setting parameters
    altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-45:0;            # (°) angle-limits for the electron beams
    E_max = 100;                   # (eV) upper limit to the energy grid
    B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith
    msis_file = find_msis_file();
    iri_file = find_iri_file();
    savedir = make_savedir("", "")
    # Define input parameters
    input_type = "LET"
    IeE_tot = 1e-3              # total energy flux (W/m²)
    E0 = 50                     # characteristic energy (eV)
    Beams = 1                   # indices of the beams with precipitation
    low_energy_tail = true      # with or without a low energy tail
    INPUT_OPTIONS = (;input_type, IeE_tot, E0, Beams, low_energy_tail);
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
    msis_file = find_msis_file();
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
    θ_lims = 180:-45:0;             # (°) angle-limits for the electron beams
    E_max = 100;                    # (eV) upper limit to the energy grid
    B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith
    t_total = 0.1;                  # (s) total simulation time
    dt = 0.01;                      # (s) time step for saving data
    n_loop = 2;                     # number of loops to run
    CFL_number = 128;
    msis_file = find_msis_file();
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
    calculate_e_transport(altitude_lims, θ_lims, E_max, B_angle_to_zenith, t_total, dt,
                          msis_file, iri_file, savedir, INPUT_OPTIONS, CFL_number; n_loop)
    @test true
end

@testitem "(TD) Does flickering run?" begin
    # Setting parameters
    altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-45:0;             # (°) angle-limits for the electron beams
    E_max = 100;                    # (eV) upper limit to the energy grid
    B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith
    t_total = 0.1;                  # (s) total simulation time
    dt = 0.01;                      # (s) time step for saving data
    n_loop = 2;                     # number of loops to run
    CFL_number = 128;
    msis_file = find_msis_file();
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
    calculate_e_transport(altitude_lims, θ_lims, E_max, B_angle_to_zenith, t_total, dt,
                          msis_file, iri_file, savedir, INPUT_OPTIONS, CFL_number; n_loop)
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
