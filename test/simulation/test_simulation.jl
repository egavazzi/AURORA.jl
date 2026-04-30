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
        @test size(sim.cache.Ie, 2) == sim.time.max_loop_length
        @test size(sim.cache.Ie_save, 2) == sim.time.max_save_steps
        @test size(sim.cache.Ie_top, 2) == length(sim.time.t)
    end
end

@testitem "Time-dependent loop slices cover resolved and save grids" begin
    altitude_lims = [100, 200]
    θ_lims = 180:-90:0
    E_max = 100
    B_angle_to_zenith = 13

    model = AuroraModel(altitude_lims, θ_lims, E_max,
                        find_msis_file(), find_iri_file(),
                        B_angle_to_zenith)
    time = RefinedTimeGrid(model, TimeDependentMode(duration = 0.1, dt = 0.01,
                                                    CFL_number = 128, n_loop = 3))

    resolved_coverage = Int[]
    save_coverage = Int[]
    for (i_slice, slice) in enumerate(time.loop_slices)
        append!(resolved_coverage, i_slice == 1 ? collect(slice.resolved_indices) : collect(slice.resolved_indices[2:end]))
        append!(save_coverage, collect(time.save_indices[slice.output_indices]))
        @test length(slice.local_save_indices) == length(slice.output_indices)
    end

    @test resolved_coverage == collect(1:length(time.t))
    @test save_coverage == collect(time.save_indices)
end

@testitem "Time-dependent run saves disjoint save-grid slices for nondivisible loops" begin
    using MAT

    mktempdir() do savedir
        altitude_lims = [100, 200]
        θ_lims = 180:-90:0
        E_max = 100
        B_angle_to_zenith = 13

        model = AuroraModel(altitude_lims, θ_lims, E_max,
                            find_msis_file(), find_iri_file(),
                            B_angle_to_zenith)
        flux = InputFlux(FlatSpectrum(1.0; E_min=50.0), SinusoidalFlickering(5.0); beams=1:2)
        sim = AuroraSimulation(model, flux, savedir;
                               mode=TimeDependentMode(duration = 0.1, dt = 0.01,
                                                      CFL_number = 128, n_loop = 3))

        run!(sim)

        files = AURORA.list_result_files(savedir)
        t_all = Float64[]
        for file in files
            data = matread(file)
            t_run = vec(data["t_run"])
            @test issorted(t_run)
            @test size(data["Ie_ztE"], 2) == length(t_run)
            append!(t_all, t_run)
        end

        @test t_all == sim.time.save_t
        @test length(unique(t_all)) == length(t_all)
        @test vec(matread(joinpath(savedir, "Ie_incoming.mat"))["t_top"]) == collect(sim.time.t)
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

@testitem "TimeDependentMode rejects too many loops" begin
    model = AuroraModel([100, 200], 180:-90:0, 100, find_msis_file(), find_iri_file(), 13)
    @test_throws ErrorException RefinedTimeGrid(model, TimeDependentMode(duration=0.04, dt=0.01,
                                                                         CFL_number=64, n_loop=10_000))
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
