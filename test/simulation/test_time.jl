@testmodule TimeTestModel begin
    using AURORA
    const msis_file = find_msis_file()
    const iri_file  = find_iri_file()
    const model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 13)
end


# ──────────────────────────────────────────────────────────────────────────────
# Unit tests — time grid properties (no sim run)
# ──────────────────────────────────────────────────────────────────────────────

@testitem "RefinedTimeGrid: basic properties" setup=[TimeTestModel] begin
    using AURORA
    flux = InputFlux(FlatSpectrum(1.0; E_min=50.0); beams=1:2)

    # duration=0.05s, dt=0.01s, n_loop = 3
    # cld(5, 3) = 2  →  loops should get 2, 2 and 1 save intervals
    sim  = AuroraSimulation(TimeTestModel.model, flux, mktempdir();
                            mode=TimeDependentMode(duration=0.05, dt=0.01,
                                                   CFL_number=128, n_loop=3))
    time = sim.time

    @test time isa RefinedTimeGrid
    @test time.n_save == 5
    @test time.n_loop == 3
    @test time.n_save_per_loop == 2
    @test time.n_t_per_loop    == 2 * time.CFL_factor + 1

    # dt_internal is exact: dt / CFL_factor
    @test time.dt_internal   ≈ time.dt / time.CFL_factor  atol=1e-15

    # Length of grids
    @test length(time.t)      == time.n_save * time.CFL_factor + 1
    @test length(time.t_save) == time.n_save + 1

    # t_save is an exact integer-stride subset of the internal grid
    @test collect(time.t[1:time.CFL_factor:end]) ≈ collect(time.t_save)  atol=1e-15

    # Internal grid starts at 0 and ends at duration
    @test time.t[begin]      ≈ 0.0        atol=1e-15
    @test time.t[end]        ≈ time.duration  atol=1e-15
    @test time.t_save[begin] ≈ 0.0        atol=1e-15
    @test time.t_save[end]   ≈ time.duration  atol=1e-15
end

@testitem "RefinedTimeGrid: loop_save_count gives correct partition" setup=[TimeTestModel] begin
    using AURORA
    flux = InputFlux(FlatSpectrum(1.0; E_min=50.0); beams=1:2)

    # Non-divisible: n_save=5, n_loop=3  →  [2, 2, 1]
    sim_nd = AuroraSimulation(TimeTestModel.model, flux, mktempdir();
                              mode=TimeDependentMode(duration=0.05, dt=0.01,
                                                     CFL_number=128, n_loop=3))
    t_nd = sim_nd.time
    @test AURORA.loop_save_count(t_nd, 1) == 2
    @test AURORA.loop_save_count(t_nd, 2) == 2
    @test AURORA.loop_save_count(t_nd, 3) == 1   # remainder

    # Sum equals total save intervals
    total = sum(AURORA.loop_save_count(t_nd, i) for i in 1:t_nd.n_loop)
    @test total == t_nd.n_save

    # Divisible: n_save=6, n_loop=3  →  [2, 2, 2]
    sim_d = AuroraSimulation(TimeTestModel.model, flux, mktempdir();
                             mode=TimeDependentMode(duration=0.06, dt=0.01,
                                                    CFL_number=128, n_loop=3))
    t_d = sim_d.time
    @test t_d.n_save == 6
    @test AURORA.loop_save_count(t_d, 1) == 2
    @test AURORA.loop_save_count(t_d, 2) == 2
    @test AURORA.loop_save_count(t_d, 3) == 2
end

@testitem "RefinedTimeGrid: loop_internal_start covers disjoint windows" setup=[TimeTestModel] begin
    using AURORA
    flux = InputFlux(FlatSpectrum(1.0; E_min=50.0); beams=1:2)

    sim  = AuroraSimulation(TimeTestModel.model, flux, mktempdir();
                            mode=TimeDependentMode(duration=0.05, dt=0.01,
                                                   CFL_number=128, n_loop=3))
    time = sim.time
    n_int_total = length(time.t)

    # Each loop's window starts where the previous one ended (shared boundary)
    for i in 1:time.n_loop - 1
        s_this = AURORA.loop_internal_start(time, i)
        n_this = AURORA.loop_internal_count(time, i)
        s_next = AURORA.loop_internal_start(time, i + 1)
        @test s_next == s_this + n_this - 1   # overlap by exactly 1 boundary point
    end

    # Last loop ends exactly at the end of the internal grid
    last_start = AURORA.loop_internal_start(time, time.n_loop)
    last_count = AURORA.loop_internal_count(time, time.n_loop)
    @test last_start + last_count - 1 == n_int_total
end


# ──────────────────────────────────────────────────────────────────────────────
# Integration tests — check simulation_data.nc time axis after run!
# ──────────────────────────────────────────────────────────────────────────────

@testitem "Time-dependent 2 loops: t_run is exact, no overlap, no drift" setup=[TimeTestModel] begin
    using NCDatasets
    mktempdir() do savedir
        flux = InputFlux(FlatSpectrum(1e-2; E_min=50.0), SmoothOnset(0.0, 0.05); beams=1:2)

        # duration=0.05s, dt=0.01s, n_save=5, n_loop=2
        sim = AuroraSimulation(TimeTestModel.model, flux, savedir;
                               mode=TimeDependentMode(duration=0.05, dt=0.01,
                                                      CFL_number=128, n_loop=2))
        run!(sim)

        NCDataset(joinpath(savedir, "simulation_data.nc"), "r") do ds
            t_all = Array(ds["time"])

            # Full time axis equals the expected uniform save grid
            expected = collect(range(0.0, 0.05; length=6))   # [0, 0.01, …, 0.05]
            @test length(t_all) == 6
            @test t_all ≈ expected  atol=1e-12
        end
    end
end

@testitem "Time-dependent 3 loops (non-divisible): t_run is exact" setup=[TimeTestModel] begin
    using NCDatasets
    mktempdir() do savedir
        flux = InputFlux(FlatSpectrum(1e-2; E_min=50.0), SmoothOnset(0.0, 0.05); beams=1:2)

        # duration=0.05s, dt=0.01s, n_save=5, n_loop=3
        sim = AuroraSimulation(TimeTestModel.model, flux, savedir;
                               mode=TimeDependentMode(duration=0.05, dt=0.01,
                                                      CFL_number=128, n_loop=3))
        run!(sim)

        NCDataset(joinpath(savedir, "simulation_data.nc"), "r") do ds
            t_all = Array(ds["time"])

            # Full time axis equals the expected uniform save grid
            expected = collect(range(0.0, 0.05; length=6))
            @test length(t_all) == 6
            @test t_all ≈ expected  atol=1e-12
        end
    end
end

@testitem "Time-dependent 3 loops (divisible): uniform loop sizes" setup=[TimeTestModel] begin
    using NCDatasets
    mktempdir() do savedir
        flux = InputFlux(FlatSpectrum(1e-2; E_min=50.0), SmoothOnset(0.0, 0.06); beams=1:2)

        # duration=0.06s, dt=0.01s, n_save=6, n_loop=3
        sim = AuroraSimulation(TimeTestModel.model, flux, savedir;
                               mode=TimeDependentMode(duration=0.06, dt=0.01,
                                                      CFL_number=128, n_loop=3))
        run!(sim)

        NCDataset(joinpath(savedir, "simulation_data.nc"), "r") do ds
            t_all = Array(ds["time"])

            # Full time axis equals the expected uniform save grid
            expected = collect(range(0.0, 0.06; length=7))
            @test length(t_all) == 7
            @test t_all ≈ expected  atol=1e-12
        end
    end
end
