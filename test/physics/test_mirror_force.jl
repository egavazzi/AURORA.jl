@testitem "DipoleField profile" begin
    field = DipoleField()
    @test field(0.0) ≈ 1.0
    @test field(600e3) < field(100e3) < 1.0
    @test DipoleField(R_earth = 1e6)(1e6) ≈ 1 / 2^3
end

@testitem "Mirror operator structure and conservation" begin
    msis_file = find_msis_file()
    iri_file = find_iri_file()
    R_earth = 6.371e6
    model = AuroraModel([100, 400], 180:-10:0, 100, msis_file, iri_file;
                        magnetic_field = DipoleField(; R_earth))

    n_z = model.altitude_grid.n
    n_beams = model.pitch_angle_grid.n_beams
    Mmirror = zeros(n_z, n_beams, n_beams)
    AURORA.update_Mmirror!(Mmirror, model)

    # Beams couple only to themselves and to their lower-μ neighbour (the upwind side:
    # the adiabatic μ-drift points toward μ = +1 everywhere in a field converging downward)
    for k in 1:n_beams, j in 1:n_beams
        if j != k && j != k - 1
            @test all(Mmirror[:, k, j] .== 0)
        end
    end
    # Upwind gains enter as negative couplings (the operator joins Mlhs with a + sign)
    for k in 2:n_beams
        @test all(Mmirror[:, k, k - 1] .< 0)
    end

    # Column sums telescope to the flux-tube convergence term κ·μ̄ⱼ: the redistribution
    # between beams cancels exactly, conserving particle number in the flux-tube sense.
    # For a vertical dipole field line κ = 3/(R_earth + h).
    h = model.altitude_grid.h
    μ_center = model.pitch_angle_grid.μ_center
    for iz in (2, n_z ÷ 2, n_z - 1)
        κ = 3 / (R_earth + h[iz])
        for j in 1:n_beams
            @test sum(Mmirror[iz, :, j]) ≈ κ * μ_center[j] rtol = 1e-3
        end
    end

    # A tilted field line stretches s by 1/cos(angle), weakening the gradient along s
    tilted = AuroraModel([100, 400], 180:-10:0, 100, msis_file, iri_file, 60;
                         magnetic_field = DipoleField(; R_earth))
    M_tilted = zeros(n_z, n_beams, n_beams)
    AURORA.update_Mmirror!(M_tilted, tilted)
    @test M_tilted ≈ cosd(60) .* Mmirror rtol = 1e-12
end

@testitem "Mirror operator off, uniform field, invalid field" begin
    uniform_field(h) = 1.0
    increasing_field(h) = exp(h / 1e6)

    msis_file = find_msis_file()
    iri_file = find_iri_file()
    model = AuroraModel([100, 400], 180:-10:0, 100, msis_file, iri_file)

    # No magnetic field profile → operator stays zero (mirroring off)
    n_z = model.altitude_grid.n
    n_beams = model.pitch_angle_grid.n_beams
    Mmirror = ones(n_z, n_beams, n_beams)
    AURORA.update_Mmirror!(Mmirror, model)
    @test all(Mmirror .== 0)

    # Uniform field → κ = 0 → zero operator (and bit-identical results to mirroring off)
    model.magnetic_field = uniform_field
    fill!(Mmirror, 1.0)
    AURORA.update_Mmirror!(Mmirror, model)
    @test all(Mmirror .== 0)

    # A field increasing with altitude breaks the upwind direction → rejected
    model.magnetic_field = increasing_field
    @test_throws ArgumentError AURORA.update_Mmirror!(Mmirror, model)

    # Bare anonymous functions are rejected at construction (not reproducible through JLD2)
    @test_throws ArgumentError AuroraModel([100, 400], 180:-10:0, 100, msis_file, iri_file;
                                           magnetic_field = h -> 1.0)
end

@testitem "Steady state with mirror force" begin
    msis_file = find_msis_file()
    iri_file = find_iri_file()
    θ_lims = 180:-90:0
    flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1)

    model_mirror = AuroraModel([100, 200], θ_lims, 100, msis_file, iri_file;
                               magnetic_field = DipoleField())
    dir_mirror = mktempdir()
    run!(AuroraSimulation(model_mirror, flux, dir_mirror; mode = SteadyStateMode()))
    Ie_mirror = load_results(dir_mirror).Ie

    model_plain = AuroraModel([100, 200], θ_lims, 100, msis_file, iri_file)
    dir_plain = mktempdir()
    run!(AuroraSimulation(model_plain, flux, dir_plain; mode = SteadyStateMode()))
    Ie_plain = load_results(dir_plain).Ie

    @test all(isfinite, Ie_mirror)
    @test all(Ie_mirror .>= 0)
    # Over 100-200 km the mirror force is a small correction, but it must enter the solve
    @test Ie_mirror != Ie_plain
    @test Ie_mirror ≈ Ie_plain rtol = 0.1
end
