@testitem "DipoleField profile" begin
    field = DipoleField()
    @test field(0.0) ≈ 1.0
    @test field(600e3) < field(100e3) < 1.0
    @test DipoleField(R_earth = 1e6)(1e6) ≈ 1 / 2^3
end

@testitem "mirror_kappa on DipoleField" begin
    R_earth = 6.371e6
    field = DipoleField(; R_earth)
    h = collect(range(100e3, 400e3, length = 50))
    r = R_earth .+ h

    # Polar (radial) field line: κ = 3/r exactly (r_eq = ∞, no latitude term)
    κ_polar = AURORA.mirror_kappa(field, h, 0)
    @test κ_polar ≈ 3 ./ r rtol = 1e-12
    @test all(κ_polar .> 0)
    @test all(diff(κ_polar) .< 0)   # decreasing with altitude

    # Tilted field line: κ stays finite, positive and decreasing with altitude
    κ_tilted = AURORA.mirror_kappa(field, h, 13)
    @test all(isfinite, κ_tilted)
    @test all(κ_tilted .> 0)
    @test all(diff(κ_tilted) .< 0)

    # Near-horizontal field line crosses the magnetic equator inside the domain → rejected
    @test_throws ArgumentError AURORA.mirror_kappa(field, h, 89)
end

@testitem "fill_mirror_operator! pure core" begin
    θ_lims = 180:-30:0
    μ_lims = cosd.(θ_lims)                 # ascending from −1 (180°) to +1 (0°)
    μ_center = (μ_lims[1:end-1] .+ μ_lims[2:end]) ./ 2
    Ω = 2π .* diff(μ_lims)
    n_beams = length(μ_center)

    n_z = 4
    κ = fill(2.0, n_z)                     # constant inverse focusing length
    Mmirror = zeros(n_z, n_beams, n_beams)
    AURORA.fill_mirror_operator!(Mmirror, μ_lims, μ_center, Ω, κ)

    # (a) Coupling structure: each beam couples only to itself and to its lower-μ neighbour
    for k in 1:n_beams, j in 1:n_beams
        if j != k && j != k - 1
            @test all(Mmirror[:, k, j] .== 0)
        end
    end
    # Upwind gains from the beam below in μ enter as non-positive couplings
    for k in 2:n_beams
        @test all(Mmirror[:, k, k - 1] .<= 0)
    end
    # No flux through the μ = ±1 edges: the lowest beam (μ = −1 edge) only self-couples, with
    # no gain from below; the top beam (μ = +1 edge) has no advective loss across its upper
    # edge — its diagonal is the pure convergence term κ·μ̄, since (1 − μ_lims[end]²) = 0.
    for j in 2:n_beams
        @test all(Mmirror[:, 1, j] .== 0)
    end
    @test all(Mmirror[iz, n_beams, n_beams] ≈ κ[iz] * μ_center[n_beams] for iz in 1:n_z)

    # (b) Flux-tube conservation: the advection telescopes so each column j sums to κ·μ̄ⱼ
    for iz in 1:n_z, j in 1:n_beams
        @test sum(Mmirror[iz, :, j]) ≈ κ[iz] * μ_center[j]
    end
end

@testitem "Mirror operator structure and conservation" begin
    msis_file = joinpath(@__DIR__, "..", "regression", "reference_results",
                         "msis_20051008-2200_70N-19E.txt")
    iri_file = joinpath(@__DIR__, "..", "regression", "reference_results",
                        "iri_20051008-2200_70N-19E.txt")
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
    # For a vertical dipole field line κ = 3/(R_earth + h), now in closed form (tight tol).
    h = model.altitude_grid.h
    μ_center = model.pitch_angle_grid.μ_center
    for iz in (2, n_z ÷ 2, n_z - 1)
        κ = 3 / (R_earth + h[iz])
        for j in 1:n_beams
            @test sum(Mmirror[iz, :, j]) ≈ κ * μ_center[j] rtol = 1e-12
        end
    end

    # A tilted field line stretches s by 1/cos(angle), but the closed-form κ also picks up the
    # latitude-dependent dipole term, so it is *more* than the pure cos-scaled radial value.
    tilted = AuroraModel([100, 400], 180:-10:0, 100, msis_file, iri_file, 60;
                         magnetic_field = DipoleField(; R_earth))
    M_tilted = zeros(n_z, n_beams, n_beams)
    AURORA.update_Mmirror!(M_tilted, tilted)

    # Expected tilted κ from the closed form (dip = 30°, anchored at the domain bottom)
    r = R_earth .+ h
    r_eq = r[1] / cosd(atand(tand(30) / 2))^2
    κ_tilted = [cosd(60) * (3 / r[iz] + 1.5 / (r_eq * (1 + 3 * (1 - r[iz] / r_eq))))
                for iz in 1:n_z]
    for iz in (2, n_z ÷ 2, n_z - 1)
        for j in 1:n_beams
            @test sum(M_tilted[iz, :, j]) ≈ κ_tilted[iz] * μ_center[j] rtol = 1e-12
        end
        # The latitude term makes it exceed the pure cos-scaled radial κ = cos(60)·3/r
        @test κ_tilted[iz] > cosd(60) * 3 / r[iz]
    end
end

@testitem "Mirror operator off, uniform field, invalid field" begin
    # Tiny functor profiles (bare functions are no longer accepted as `magnetic_field`)
    struct UniformField <: AbstractFieldProfile end
    (::UniformField)(h) = 1.0
    struct IncreasingField <: AbstractFieldProfile end
    (::IncreasingField)(h) = exp(h / 1e6)

    msis_file = joinpath(@__DIR__, "..", "regression", "reference_results",
                         "msis_20051008-2200_70N-19E.txt")
    iri_file = joinpath(@__DIR__, "..", "regression", "reference_results",
                        "iri_20051008-2200_70N-19E.txt")
    model = AuroraModel([100, 400], 180:-10:0, 100, msis_file, iri_file)

    # No magnetic field profile → operator stays zero (mirroring off)
    n_z = model.altitude_grid.n
    n_beams = model.pitch_angle_grid.n_beams
    Mmirror = ones(n_z, n_beams, n_beams)
    AURORA.update_Mmirror!(Mmirror, model)
    @test all(Mmirror .== 0)

    # Uniform field → κ = 0 (default finite-difference path) → zero operator
    model.magnetic_field = UniformField()
    fill!(Mmirror, 1.0)
    AURORA.update_Mmirror!(Mmirror, model)
    @test all(Mmirror .== 0)

    # A field increasing with altitude breaks the upwind direction → rejected
    model.magnetic_field = IncreasingField()
    @test_throws ArgumentError AURORA.update_Mmirror!(Mmirror, model)

    # Bare anonymous functions are rejected at construction (not reproducible through JLD2)
    @test_throws ArgumentError AuroraModel([100, 400], 180:-10:0, 100, msis_file, iri_file;
                                           magnetic_field = h -> 1.0)
end

@testitem "Steady state with mirror force" begin
    msis_file = joinpath(@__DIR__, "..", "regression", "reference_results",
                         "msis_20051008-2200_70N-19E.txt")
    iri_file = joinpath(@__DIR__, "..", "regression", "reference_results",
                        "iri_20051008-2200_70N-19E.txt")
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
