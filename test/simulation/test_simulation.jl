@testitem "AuroraSimulation initialize! populates cache" begin
    mktempdir() do savedir
        altitude_lims = [100, 400]
        θ_lims = 180:-45:0
        E_max = 100
        B_angle_to_zenith = 13
        t_total = 0.1
        dt = 0.01
        CFL_number = 128
        msis_file = find_msis_file()
        iri_file = find_iri_file()

        model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)
        flux = InputFlux(FlatSpectrum(1.0; E_min=50.0), SmoothOnset(0.0, 0.05);
                         beams=1:2, z_source=500.0)

        sim = AuroraSimulation(model, flux, savedir;
                               solver=TimeDependentSolver(t_total, dt; CFL_number, n_loop=2))

        @test sim.cache === nothing
        @test sim.time isa RefinedTimeGrid
        @test sim.time.dt_resolved <= sim.time.dt_requested

        initialize!(sim)

        @test sim.cache !== nothing
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
        sim = AuroraSimulation(model, flux, savedir; solver=SteadyStateSolver())

        @test sim.cache === nothing

        run!(sim)

        @test sim.cache !== nothing
    end
end
