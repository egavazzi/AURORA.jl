@testitem "AuroraSimulation initialize! populates cache" begin
    mktempdir() do savedir
        altitude_lims = [100, 400]
        θ_lims = 180:-45:0
        E_max = 100
        B_angle_to_zenith = 13

        msis_file = find_msis_file()
        iri_file = find_iri_file()

        model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)
        flux = InputFlux(FlatSpectrum(1.0; E_min=50.0), SmoothOnset(0.0, 0.05);
                         beams=1:2, z_source=500.0)

        sim = AuroraSimulation(model, flux, savedir;
                               mode=TimeDependentMode(duration = 0.1, dt = 0.01,
                                                      CFL_number = 128, n_loop = 2))

        @test !sim.cache_initialized
        @test sim.time isa RefinedTimeGrid
        @test sim.time.dt_internal <= sim.time.dt

        initialize!(sim)

        @test sim.cache_initialized
        @test sim.model.initialized
        n_species = 3
        @test sim.cache.degradation.secondary_e_flux isa NTuple{n_species, Matrix{Float64}}
        @test sim.cache.degradation.primary_e_spectrum isa NTuple{n_species, Vector{Float64}}
        @test all(sp.cascading_data isa AURORA.SpeciesCascadingCache for sp in sim.model.species)
        @test all(!isempty(sp.cascading_data.E_edges) for sp in sim.model.species)
        @test size(sim.cache.Ie, 2) == sim.time.n_t_per_loop
        @test size(sim.cache.Ie_top, 2) == length(sim.time.t)
    end
end

@testitem "AuroraSimulation run! auto-initializes" begin
    mktempdir() do savedir
        altitude_lims = [100, 400]
        θ_lims = 180:-45:0
        E_max = 100
        B_angle_to_zenith = 13
        msis_file = find_msis_file()
        iri_file = find_iri_file()

        model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)
        flux = InputFlux(FlatSpectrum(1.0; E_min=50.0); beams=1:2)
        sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode())

        @test !sim.cache_initialized

        run!(sim)

        @test sim.cache_initialized
    end
end

@testitem "Multi-step SS: saved t_run matches time grid" begin
    using MAT
    mktempdir() do savedir
        altitude_lims = [100, 200]
        θ_lims = 180:-90:0
        E_max = 100
        B_angle_to_zenith = 13

        msis_file = find_msis_file()
        iri_file = find_iri_file()

        model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)
        flux = InputFlux(FlatSpectrum(1.0; E_min=50.0), SinusoidalFlickering(5.0); beams=1:2)
        sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode(duration = 0.04, dt = 0.01))

        run!(sim)

        data = matread(joinpath(savedir, "IeFlickering-01.mat"))
        t_run = vec(data["t_run"])
        expected_t = collect(sim.time.t)

        # t_run must span the full time axis, not be a scalar 1
        @test length(t_run) == length(expected_t)
        @test t_run ≈ expected_t

        # Ie_ztE time dimension must also match
        @test size(data["Ie_ztE"], 2) == length(expected_t)
    end
end

@testitem "SteadyStateMode() → SingleStepConfig" begin
    mktempdir() do savedir
        msis_file = find_msis_file()
        iri_file  = find_iri_file()
        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 13)
        flux  = InputFlux(FlatSpectrum(1.0; E_min=50.0); beams=1:2)

        sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode())

        @test sim.time isa SingleStepConfig
        @test sim.time.n_steps == 1
        @test sim.time.t == 1:1
    end
end

@testitem "SteadyStateMode divisibility check" begin
    # duration not an integer multiple of dt → error
    @test_throws ErrorException SteadyStateMode(duration=0.05, dt=0.02)

    # exact multiple → succeeds
    mode = SteadyStateMode(duration=0.04, dt=0.01)
    @test mode.duration == 0.04
    @test mode.dt == 0.01
end

@testitem "TimeDependentMode divisibility check" begin
    # duration not an integer multiple of dt → error
    @test_throws ErrorException TimeDependentMode(duration=0.05, dt=0.02, CFL_number=64)

    # exact multiple → succeeds
    mode = TimeDependentMode(duration=0.04, dt=0.01, CFL_number=64)
    @test mode.duration == 0.04
    @test mode.dt == 0.01
end

@testitem "Mode aliases construct canonical types" begin
    steady_state = SteadyState()
    multi_step = SteadyState(duration=0.04, dt=0.01)
    time_dependent = TimeDependent(duration=0.04, dt=0.01, CFL_number=64)

    @test steady_state isa SteadyStateMode
    @test multi_step isa SteadyStateMode
    @test time_dependent isa TimeDependentMode

    @test isnothing(steady_state.duration)
    @test multi_step.duration == 0.04
    @test multi_step.dt == 0.01
    @test time_dependent.duration == 0.04
    @test time_dependent.dt == 0.01
end

@testitem "NeutralSpecies density_profile types" begin
    msis_file = find_msis_file()
    iri_file  = find_iri_file()
    model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)

    # Default model: all species backed by MSISDensity; density is empty before initialize!
    for sp in model.species
        @test sp.density_profile isa MSISDensity
    end
    @test isempty(model.species[1].density)

    # After initialize! density is populated
    initialize!(model)
    ag = model.altitude_grid
    @test !isempty(model.species[1].density)

    # VectorDensity round-trip: PCHIP through exact sample points → same density as MSISDensity
    raw_n2 = AURORA.load_msis_density(msis_file, :N2, ag.h)
    vd = VectorDensity(ag.h, raw_n2)
    model_vd = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
    model_vd.species[1].density_profile = vd
    initialize!(model_vd)
    @test model_vd.species[1].density_profile isa VectorDensity
    @test model_vd.species[1].density ≈ model.species[1].density rtol=1e-6

    # Plain callable is also accepted as density_profile; density is empty until initialize!
    flat_profile = h -> fill(1e15, length(h))
    sp_fn = AURORA.N2Species(flat_profile)
    @test sp_fn.density_profile === flat_profile
    @test isempty(sp_fn.density)
end

@testitem "AuroraModel species support Symbol indexing" begin
    msis_file = find_msis_file()
    iri_file  = find_iri_file()
    model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)

    @test model.species[:N2] === model.species[1]
    @test model.species[:O2] === model.species[2]
    @test model.species[:O] === model.species[3]
    @test_throws KeyError model.species[:NO]
end

@testitem "Species Symbol indexing rejects duplicate names" begin
    msis_file = find_msis_file()
    iri_file  = find_iri_file()
    model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0;
                        species = (N2Species(msis_file), N2Species(msis_file)))

    @test_throws ArgumentError model.species[:N2]
end

@testitem "AuroraModel is uninitialized before initialize!" begin
    msis_file = find_msis_file()
    iri_file  = find_iri_file()
    model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)

    @test !model.initialized
    @test model.scattering isa AURORA.ScatteringData
    @test isempty(model.scattering.θ_scatter)
    @test isempty(model.species[1].density)
end

@testitem "initialize!(model) interception window" begin
    mktempdir() do savedir
        msis_file = find_msis_file()
        iri_file  = find_iri_file()
        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)

        flat_n2 = h -> fill(1e18, length(h))
        model.species[1].density_profile = flat_n2

        flux = InputFlux(FlatSpectrum(1.0; E_min=50.0); beams=1:2)
        sim  = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode())

        @test !sim.model.initialized
        run!(sim)
        @test sim.model.initialized

        n2_density = sim.model.species[1].density
        @test !isempty(n2_density)
        @test n2_density[1] ≈ 1e18   # bottom of the grid, unaffected by boundary taper
    end
end

@testitem "AuroraModel with 2 species: run! succeeds" begin
    mktempdir() do savedir
        msis_file = find_msis_file()
        iri_file  = find_iri_file()

        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0;
                            species = (O2Species(msis_file), OSpecies(msis_file)))
        flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim  = AuroraSimulation(model, flux, savedir; mode = SteadyStateMode())
        run!(sim)

        @test sim.model.initialized
        @test length(sim.model.species) == 2
        @test sim.cache.degradation.secondary_e_flux isa NTuple{2, Matrix{Float64}}
    end
end

@testitem "Custom 4th species with pre-populated cross sections: run! succeeds" begin
    mktempdir() do savedir
        msis_file = find_msis_file()
        iri_file  = find_iri_file()

        custom_law  = (E_s, E_p) -> 1.0 / (11.4^2 + E_s^2)
        custom_spec = AURORA.CascadingSpec("CustomGas", [15.581, 16.73, 18.75], custom_law)
        custom_sp   = AURORA.NeutralSpecies(:CustomGas, h -> fill(1e18, length(h));
                                            cascading_spec      = custom_spec,
                                            phase_fcn_generator = AURORA.phase_fcn_N2)

        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0;
                            species = (N2Species(msis_file), O2Species(msis_file),
                                       OSpecies(msis_file), custom_sp))

        # Interception window: pre-populate cross sections and excitation levels before
        # initialize!(model) runs (name-based auto-lookup would fail for :CustomGas)
        # We just reuse the N2 cross sections and levels here
        eg = model.energy_grid
        model.species[end].cross_sections    = AURORA.get_cross_section("N2", eg.E_centers)
        model.species[end].excitation_levels = AURORA.load_excitation_threshold_for("N2")

        flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim  = AuroraSimulation(model, flux, savedir; mode = SteadyStateMode())
        run!(sim)

        @test sim.model.initialized
        @test length(sim.model.species) == 4
        @test sim.cache.degradation.secondary_e_flux isa NTuple{4, Matrix{Float64}}
        @test sim.model.species[end].name == :CustomGas
        @test !isempty(sim.model.species[end].density)
    end
end

@testitem "Custom phase function via interception window: run! succeeds" begin
    mktempdir() do savedir
        msis_file = find_msis_file()
        iri_file  = find_iri_file()

        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)

        custom_generator = (θ, E) -> AURORA.phase_fcn_N2(θ, E)
        model.species[1].phase_fcn_generator = custom_generator

        flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim  = AuroraSimulation(model, flux, savedir; mode = SteadyStateMode())
        run!(sim)

        @test sim.model.initialized
        @test model.species[1].phase_fcn_generator === custom_generator
    end
end

@testitem "Altitude grid swap: initialize!(model) rebuilds s_field and ionosphere" begin
    msis_file = find_msis_file()
    iri_file  = find_iri_file()

    model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
    initialize!(model)
    old_n_z = model.altitude_grid.n
    @test length(model.s_field)               == old_n_z
    @test length(model.ionosphere.Tn)         == old_n_z
    @test length(model.species[1].density)    == old_n_z

    model.altitude_grid = AltitudeGrid(100, 300)
    initialize!(model)

    new_n_z = model.altitude_grid.n
    @test new_n_z > old_n_z
    @test length(model.s_field)               == new_n_z
    @test length(model.ionosphere.Tn)         == new_n_z
    @test length(model.species[1].density)    == new_n_z
end

@testitem "Altitude grid swap: run! after initialize!(model) succeeds" begin
    mktempdir() do savedir
        msis_file = find_msis_file()
        iri_file  = find_iri_file()

        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
        flux  = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim   = AuroraSimulation(model, flux, savedir; mode = SteadyStateMode())
        run!(sim)

        model.altitude_grid = AltitudeGrid(100, 300)
        initialize!(model)   # recomputes s_field, ionosphere, species
        initialize!(sim)     # rebuilds simulation cache for new grid dimensions
        run!(sim)

        @test sim.model.initialized
        @test length(sim.model.s_field) == model.altitude_grid.n
    end
end

@testitem "Reassigning a grid invalidates the model" begin
    msis_file = find_msis_file()
    iri_file  = find_iri_file()

    model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
    initialize!(model)
    @test model.initialized

    # Each geometry reassignment must flip `initialized` back to false.
    model.altitude_grid = AltitudeGrid(100, 300)
    @test !model.initialized

    initialize!(model)
    model.energy_grid = EnergyGrid(200)
    @test !model.initialized

    initialize!(model)
    model.B_angle_to_zenith = 20
    @test !model.initialized
end

@testitem "Grid change then bare run! auto-reinitializes (no manual init)" begin
    mktempdir() do savedir
        msis_file = find_msis_file()
        iri_file  = find_iri_file()

        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
        flux  = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim   = AuroraSimulation(model, flux, savedir; mode = SteadyStateMode())
        run!(sim)

        # Change the grid and call run! directly — no initialize!(model)/initialize!(sim).
        model.altitude_grid = AltitudeGrid(100, 300)
        run!(sim)

        @test sim.model.initialized
        @test length(sim.model.s_field)           == model.altitude_grid.n
        @test size(sim.cache.Ie, 1) ÷ length(model.pitch_angle_grid.μ_center) == model.altitude_grid.n
    end
end

@testitem "Energy grid change rebuilds sim.time and cache (TimeDependent)" begin
    mktempdir() do savedir
        msis_file = find_msis_file()
        iri_file  = find_iri_file()

        model = AuroraModel([100, 200], 180:-90:0, 80, msis_file, iri_file, 0)
        flux  = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim   = AuroraSimulation(model, flux, savedir;
                                 mode = TimeDependentMode(duration=0.02, dt=0.01,
                                                          CFL_number=128, n_loop=1))
        run!(sim)
        @test size(sim.cache.Ie, 3) == model.energy_grid.n

        # Larger energy grid → more energy bins AND a different CFL-refined time grid.
        model.energy_grid = EnergyGrid(200)
        run!(sim)

        @test sim.model.initialized
        @test size(sim.cache.Ie, 3) == model.energy_grid.n
        @test sim.time isa AURORA.RefinedTimeGrid
    end
end
