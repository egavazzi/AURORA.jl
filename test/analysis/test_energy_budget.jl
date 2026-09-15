@testitem "steady-state energy budget identities" setup=[SharedSimResults] begin
    using AURORA
    using NCDatasets

    dir = SharedSimResults.ss_dir
    budget = energy_budget(dir; verbose = false)

    @test budget isa EnergyBudget
    @test budget.interval === nothing          # a snapshot, so fluxes rather than energies
    @test budget.input > 0
    @test budget.escape >= 0
    @test budget.inelastic >= 0
    @test budget.heating >= 0
    @test budget.ionpairs > 0
    @test budget.net == budget.input - budget.escape
    @test budget.inelastic ≈ budget.ionization + budget.excitation
    @test budget.inelastic ≈ sum(last, budget.inelastic_by_species)
    @test budget.residual ≈ budget.input - budget.inelastic - budget.heating - budget.escape
    @test budget.residual_fraction ≈ budget.residual / budget.input
    @test budget.albedo ≈ budget.escape / budget.input
    @test [name for (name, _) in budget.inelastic_by_species] == ["N2", "O2", "O"]

    # The input matches the spectrum's normalization, converted from W/m² to eV m⁻² s⁻¹.
    @test budget.input ≈ 1e-2 / AURORA.qₑ

    # The top-altitude fluxes are the same quantity that make_current_file writes.
    make_current_file(dir)
    NCDataset(joinpath(dir, "analysis", "currents.nc"), "r") do ds
        top = argmax(ds["altitude"][:])
        @test budget.input ≈ ds["IeE_down"][top, end]
        @test budget.escape ≈ ds["IeE_up"][top, end]
    end
end

@testitem "energy budget summary" setup=[SharedSimResults] begin
    using AURORA

    budget = energy_budget(SharedSimResults.ss_dir; verbose = false)
    text = sprint(show, MIME"text/plain"(), budget)
    lines = split(text, '\n')

    @test occursin("EnergyBudget — steady state, ∫ along the field line", lines[1])
    @test occursin("eV m⁻² s⁻¹", lines[2])
    @test occursin("% of input", lines[3])
    for name in ("ionization", "excitation", "thermal heating", "escape (↑ top)", "residual")
        @test any(line -> occursin(name, line), lines)
    end
    @test occursin("inelastic by species (% of inelastic): N2 ", text)
    @test occursin("net energy per ion pair", text)
    @test occursin("(deposited + backscattered)/input", text)

    # The five percentage rows are shares of the input, so they sum to 100%.
    shares = (budget.ionization, budget.excitation, budget.heating, budget.escape,
              budget.residual)
    @test sum(shares) ≈ budget.input

    # A small but nonzero term keeps two significant digits instead of reading as "0.0".
    @test AURORA.percent_string(2e-4, 1.0) == "0.02"
    @test AURORA.percent_string(0.5, 1.0) == "50.0"
    @test AURORA.percent_string(0.0, 1.0) == "0.0"
end

@testitem "energy budget from a simulation in memory" begin
    using AURORA

    altitude_lims = [100, 200]
    θ_lims = 180:-90:0
    E_max = 100
    B_angle_to_zenith = 13
    model = AuroraModel(altitude_lims, θ_lims, E_max, find_msis_file(; verbose = false),
                        find_iri_file(; verbose = false), B_angle_to_zenith)
    savedir = mktempdir()
    sim = AuroraSimulation(model, InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2),
                           savedir; mode = SteadyStateMode())
    run!(sim; verbose = false)

    from_memory = energy_budget(sim; verbose = false)
    from_disk = energy_budget(savedir; verbose = false)
    for name in AURORA.ENERGY_BUDGET_SCALAR_FIELDS
        @test getfield(from_memory, name) ≈ getfield(from_disk, name)
    end
    @test from_memory.inelastic_by_species == from_disk.inelastic_by_species

    @test_throws "tidx = 2 out of range 1:1" energy_budget(sim; tidx = 2, verbose = false)
    @test_throws "tidx = 0 out of range 1:1" energy_budget(savedir; tidx = 0, verbose = false)
    @test_throws "pass either tidx" energy_budget(savedir; tidx = 1, trange = 1:1,
                                                  verbose = false)

    # TOML round trip, from the simulation and from its directory.
    saved = make_energy_budget_file(sim; verbose = false)
    @test isfile(joinpath(savedir, "analysis", "energy_budget.toml"))
    loaded = load_energy_budget(savedir)
    for name in fieldnames(EnergyBudget)
        @test getfield(loaded, name) == getfield(saved, name)
    end
    @test load_energy_budget(sim).input == saved.input
    @test make_energy_budget_file(savedir; verbose = false).input == saved.input
end

@testitem "time-dependent energy budget" setup=[SharedSimResults] begin
    using AURORA

    dir = SharedSimResults.td_dir

    # A single slice of a time-dependent run does not close, and says so.
    redirect_stdout(devnull) do
        @test_logs (:warn, r"single-time snapshot") energy_budget(dir)
    end
    @test_logs energy_budget(dir; verbose = false)   # quiet: no summary, no warning
    snapshot = energy_budget(dir; tidx = 1, verbose = false)
    @test snapshot.input == 0             # the pulse has not reached the top yet
    @test isnan(snapshot.albedo)

    integrated = energy_budget(dir; trange = :, verbose = false)
    @test integrated.interval !== nothing
    @test integrated.interval[2] > integrated.interval[1]
    @test integrated.input > 0
    @test integrated.net == integrated.input - integrated.escape
    @test integrated.inelastic ≈ integrated.ionization + integrated.excitation
    @test integrated.inelastic ≈ sum(last, integrated.inelastic_by_species)
    @test integrated.residual ≈
          integrated.input - integrated.inelastic - integrated.heating - integrated.escape
    @test integrated.residual_fraction ≈ integrated.residual / integrated.input
    @test sprint(show, MIME"text/plain"(), integrated) |> text -> occursin("∫ over t =", text)

    # The chunk size only sets peak memory; one slice per chunk gives the same answer.
    streamed = energy_budget(dir; trange = :, max_bytes = 1, verbose = false)
    for name in AURORA.ENERGY_BUDGET_SCALAR_FIELDS
        @test getfield(streamed, name) == getfield(integrated, name)
    end

    # A sub-range integrates less energy over a shorter span.
    partial = energy_budget(dir; trange = 2:5, verbose = false)
    @test partial.interval[2] - partial.interval[1] <
          integrated.interval[2] - integrated.interval[1]
    @test partial.input < integrated.input

    # The integrated budget round-trips through TOML, interval included.
    make_energy_budget_file(dir; trange = :, verbose = false)
    loaded = load_energy_budget(dir)
    for name in fieldnames(EnergyBudget)
        @test getfield(loaded, name) == getfield(integrated, name)
    end
end

@testitem "energy budget input validation" setup=[SharedSimResults] begin
    using AURORA

    empty_dir = mktempdir()
    @test_throws "no physics_state.jld2 found" energy_budget(empty_dir; verbose = false)
    @test_throws "no energy_budget.toml found" load_energy_budget(empty_dir)

    dir = SharedSimResults.td_dir
    @test_throws "trange must be a Colon" energy_budget(dir; trange = [1, 3], verbose = false)
    @test_throws "trange must be a Colon" energy_budget(dir; trange = 0:3, verbose = false)
    @test_throws "need ≥ 2 time slices" energy_budget(dir; trange = 2:2, verbose = false)
end
