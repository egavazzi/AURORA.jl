@testitem "steady-state energy budget identities" setup=[SharedSimResults] begin
    using AURORA
    using NCDatasets

    dir = SharedSimResults.ss_dir
    budget = energy_budget(dir; verbose = false)

    @test budget isa EnergyBudget
    @test budget.interval === nothing          # a snapshot, so fluxes rather than energies
    @test budget.input > 0
    @test budget.escape >= 0
    @test budget.bottom_escape >= 0
    @test budget.inelastic >= 0
    @test budget.heating >= 0
    @test budget.ionpairs > 0
    @test budget.net == budget.input - budget.escape - budget.bottom_escape
    @test budget.inelastic ≈ budget.ionization + budget.excitation
    @test budget.inelastic ≈ sum(last, budget.inelastic_by_species)
    @test budget.subfloor > 0
    @test budget.residual ≈ budget.input - budget.inelastic - budget.heating -
                            budget.subfloor - budget.escape - budget.bottom_escape
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
    for name in ("ionization", "excitation", "thermal heating", "sub-floor thermalisation",
                 "escape (↑ top)", "escape (↓ bottom)", "residual")
        @test any(line -> occursin(name, line), lines)
    end
    @test occursin("inelastic by species (% of inelastic): N2 ", text)
    @test occursin("net energy per ion pair", text)
    @test occursin("(deposited + escaped)/input", text)

    # The seven percentage rows are shares of the input, so they sum to 100%.
    shares = (budget.ionization, budget.excitation, budget.heating, budget.subfloor,
              budget.escape, budget.bottom_escape, budget.residual)
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
    @test_throws "pass at most one of" energy_budget(savedir; tidx = 1, trange = 1:1,
                                                     verbose = false)

    # TOML round trip, from the simulation and from its directory.
    saved = make_energy_budget_file(sim; verbose = false)
    @test isfile(joinpath(savedir, "analysis", "energy_budget.toml"))
    loaded = load_energy_budget(savedir)
    for name in fieldnames(EnergyBudget)
        @test getfield(loaded, name) == getfield(saved, name)
    end
    @test loaded.subfloor == saved.subfloor
    @test load_energy_budget(sim).input == saved.input
    @test make_energy_budget_file(savedir; verbose = false).input == saved.input

    # A file without one of the budget's terms is refused rather than read with a zero.
    path = joinpath(savedir, "analysis", "energy_budget.toml")
    data = AURORA.TOML.parsefile(path)
    delete!(data["values"], "subfloor")
    open(io -> AURORA.TOML.print(io, data), path, "w")
    @test_throws "has no value for subfloor" load_energy_budget(savedir)
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
    @test integrated.net == integrated.input - integrated.escape - integrated.bottom_escape
    @test integrated.inelastic ≈ integrated.ionization + integrated.excitation
    @test integrated.inelastic ≈ sum(last, integrated.inelastic_by_species)
    @test integrated.residual ≈ integrated.input - integrated.inelastic -
                                integrated.heating - integrated.subfloor -
                                integrated.escape - integrated.bottom_escape
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

    # A (t0, t1) tuple selects the slices inside that closed time interval.
    t = load_coordinates(dir).t
    by_time = energy_budget(dir; trange = (t[2], t[5]), verbose = false)
    by_index = energy_budget(dir; trange = 2:5, verbose = false)
    for name in fieldnames(EnergyBudget)
        @test getfield(by_time, name) == getfield(by_index, name)
    end
    @test by_time.interval == (t[2], t[5])

    @test_throws "at least 2 are needed" energy_budget(dir; trange = (t[2], t[2]),
                                                       verbose = false)
    @test_throws "at least 2 are needed" energy_budget(dir; trange = (-2.0, -1.0),
                                                       verbose = false)
    @test_throws "runs backwards" energy_budget(dir; trange = (t[5], t[2]), verbose = false)

    # `t` picks the saved slice nearest that time.
    by_t = energy_budget(dir; t = t[4], verbose = false)
    at_index = energy_budget(dir; tidx = 4, verbose = false)
    for name in fieldnames(EnergyBudget)
        @test getfield(by_t, name) == getfield(at_index, name)
    end
    nearer_to_4 = energy_budget(dir; t = t[4] + 0.3 * (t[5] - t[4]), verbose = false)
    @test nearer_to_4.input == at_index.input
    nearer_to_5 = energy_budget(dir; t = t[4] + 0.7 * (t[5] - t[4]), verbose = false)
    @test nearer_to_5.input == energy_budget(dir; tidx = 5, verbose = false).input

    @test_throws "outside the saved time span" energy_budget(dir; t = last(t) + 1,
                                                             verbose = false)
    @test_throws "outside the saved time span" energy_budget(dir; t = first(t) - 1,
                                                             verbose = false)
    @test_throws "pass at most one of" energy_budget(dir; t = t[2], tidx = 2,
                                                     verbose = false)

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
    @test_throws "out of bounds for the" energy_budget(dir; trange = 0:3, verbose = false)
    @test_throws "need ≥ 2 time slices" energy_budget(dir; trange = 2:2, verbose = false)
end

@testitem "sub-floor factor of a single non-ionizing channel" begin
    using AURORA

    # Three 1 eV bins above a 1 eV floor; one channel of loss L = 1.5 eV.
    E_edges = [1.0, 2.0, 3.0, 4.0]
    energy_grid = (; E_edges, E_centers = [1.5, 2.5, 3.5], ΔE = [1.0, 1.0, 1.0])
    L = 1.5
    species = (; name = :X, cross_sections = [0.0 0.0 0.0; 0.0 1e-20 1e-20],
               excitation_levels = [0.0 0.0; L 0.0], density = fill(1e17, 3),
               cascading_data = nothing)
    F = AURORA.subfloor_energy_factors(species, energy_grid)
    @test size(F) == size(species.cross_sections)
    @test all(iszero, F[1, :])                    # elastic
    # From bin 2 the electrons arrive uniformly over [0.5, 1.5] eV. The half below the 1 eV
    # floor thermalises with a mean energy of 0.75 eV.
    @test F[2, 2] ≈ 0.5 * 0.75
    # From bin 1 they arrive over [-0.5, 0.5] eV, all below the floor; the arrivals below
    # zero energy carry none.
    @test F[2, 1] ≈ 0.5 * 0.25
    # From bin 3 they arrive over [1.5, 2.5] eV, on the grid.
    @test F[2, 3] == 0

    # Coulomb: Le·Ie/ΔE[1] electrons leave bin 1, crossing the floor at E_edges[1].
    @test AURORA.coulomb_subfloor_factors(energy_grid) == [1.0, 0.0, 0.0]

    # In the budget, flux in bin 2 only: subfloor / inelastic is F[2, 2] / L. Zero electron
    # density switches the Coulomb loss off.
    z = [100e3, 101e3, 102e3]
    model = (; energy_grid, altitude_grid = (; h = z), s_field = z,
             pitch_angle_grid = (; μ_center = [-0.5, 0.5]),
             ionosphere = (; ne = zeros(3), Te = fill(1000.0, 3)), species = [species])
    Ie = zeros(3, 2, 3)
    Ie[:, :, 2] .= 1e10
    budget = AURORA.energy_budget_snapshot(model, Ie; verbose = false)
    @test budget.heating == 0
    @test budget.inelastic > 0
    @test budget.subfloor ≈ budget.inelastic * F[2, 2] / L
end

@testitem "sub-floor term on the regression model" begin
    using AURORA
    const AU = AURORA

    # The model and input of the steady-state regression test.
    input_data = joinpath(@__DIR__, "..", "regression", "input_data")
    model = AuroraModel([100, 400], 180:-30:0, 500,
                        joinpath(input_data, "msis_20051008-2200_70N-19E.txt"),
                        joinpath(input_data, "iri_20051008-2200_70N-19E.txt"), 13)
    savedir = mktempdir()
    sim = AuroraSimulation(model, InputFlux(FlatSpectrum(1e-2; E_min = 400.0); beams = 1:2),
                           savedir; mode = SteadyStateMode())
    run!(sim; verbose = false)
    budget = energy_budget(sim; verbose = false)

    # Reference: the energy below the floor of every collision's outgoing electrons, per
    # energy bin, from the arrival range (non-ionizing) or the secondary law (ionizing), and
    # the Coulomb exits from bin 1 at the floor energy.
    const eg = model.energy_grid
    const E_edges, ΔE = eg.E_edges, eg.ΔE
    const fl = E_edges[1]
    const n_E = length(ΔE)

    function gl(f, a, b; X = AU.GL16_X, W = AU.GL16_W)
        b > a || return 0.0
        m, hw = (a + b) / 2, (b - a) / 2
        return hw * sum(w * f(m + hw * x) for (x, w) in zip(X, W))
    end
    function pieces(lo, hi, bps)
        pts = sort(unique([lo; filter(p -> lo < p < hi, bps); hi]))
        return [(pts[k], pts[k + 1]) for k in 1:length(pts) - 1]
    end
    # Sub-floor energy of the ionization events of primary bin [Emin, Emax], unnormalized.
    function single_ref(law, I, Emin, Emax)
        lo = max(Emin, I)
        total = 0.0
        Emax > lo || return total
        for (a, b) in pieces(lo, Emax, [I + fl, I + 2fl])
            total += gl(a, b) do Ep
                W = Ep - I
                gl(Es -> Es * law(Es, Ep), 0.0, min(fl, W / 2)) +
                gl(Es -> (W - Es) * law(Es, Ep), max(0.0, W - fl), W / 2)
            end
        end
        return total
    end
    function double_ref(law, I, Emin, Emax)
        lo = max(Emin, I)
        total = 0.0
        Emax > lo || return total
        for (a, b) in pieces(lo, Emax, [I + fl, I + 2fl, I + 3fl])
            total += gl(a, b; X = AU.GL8_X, W = AU.GL8_W) do Ep
                W = Ep - I
                W > 0 || return 0.0
                cdf = AU.build_secondary_cdf(E_edges, W / 2, Ep, law)
                2 * gl(Es -> Es * AU.double_secondary_density((I, cdf), Es),
                       0.0, min(fl, W / 2)) +
                gl(Ed -> Ed * AU.double_primary_density_cdf((I, cdf), Ed), W / 3, min(W, fl))
            end
        end
        return total
    end
    function reference(model, Ie_omni)
        collisional = 0.0
        for sp in model.species
            lv = sp.excitation_levels
            law = AU.runtime_law(sp.cascading_data.spec.secondary_law)
            for l in axes(lv, 1)[2:end], iE in 1:n_E
                L, n_sec, σ = lv[l, 1], Int(lv[l, 2]), sp.cross_sections[l, iE]
                σ == 0 && continue
                if n_sec == 0
                    llo, lup = E_edges[iE] - L, min(E_edges[iE + 1] - L, E_edges[iE])
                    hi, lo = min(lup, fl), max(llo, 0.0)
                    per_collision = hi > lo ? (hi^2 - lo^2) / (2ΔE[iE]) : 0.0
                else
                    raw = n_sec == 1 ? single_ref(law, L, E_edges[iE], E_edges[iE + 1]) :
                                       double_ref(law, L, E_edges[iE], E_edges[iE + 1])
                    per_collision = raw / AU.event_count(sp.cascading_data, iE, L)
                end
                collisional += σ * sum(sp.density .* Ie_omni[:, iE]) * per_collision
            end
        end
        Le = AU.loss_to_thermal_electrons(eg.E_centers[1], model.ionosphere.ne,
                                          model.ionosphere.Te)
        coulomb = sum(Le .* Ie_omni[:, 1]) / ΔE[1] * fl
        return collisional, coulomb
    end

    μ = model.pitch_angle_grid.μ_center
    n_z, n_μ = length(model.altitude_grid.h), length(μ)
    w = AU.column_weights(model.s_field)
    Ie = reshape(sim.workspace.Ie[:, 1, :], n_z, n_μ, n_E)
    Ie_omni = [sum((μ[iμ] < 0 ? w.down[iz] : w.up[iz]) * Ie[iz, iμ, iE] for iμ in 1:n_μ)
               for iz in 1:n_z, iE in 1:n_E]
    collisional, coulomb = Base.invokelatest(reference, model, Ie_omni)

    @test budget.subfloor ≈ collisional + coulomb rtol = 1e-6
    # The same terms as eV m⁻² s⁻¹ for this run, to four digits.
    @test budget.input ≈ 6.2415e16 rtol = 1e-4
    @test collisional ≈ 2.0222e15 rtol = 1e-4
    @test coulomb ≈ 1.4059e14 rtol = 1e-4
    @test budget.subfloor / budget.input ≈ 0.03465 rtol = 1e-3

    # What is left is the scheme's numerical non-conservation, 0.86% of the input.
    @test budget.residual_fraction ≈ 0.0086 rtol = 0.02
end
