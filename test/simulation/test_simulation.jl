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
        @test sim.time.dt_resolved <= sim.time.dt_requested

        initialize!(sim)

        @test sim.cache_initialized
        @test sim.cache.cascading isa AURORA.CascadingCache
        n_species = 3
        @test sim.cache.degradation.secondary_e_flux isa NTuple{n_species, Matrix{Float64}}
        @test sim.cache.degradation.primary_e_spectrum isa NTuple{n_species, Vector{Float64}}
        @test all(!isempty(species_cache.E_edges) for species_cache in sim.cache.cascading)
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
