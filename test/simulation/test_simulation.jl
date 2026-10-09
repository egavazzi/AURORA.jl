@testitem "AuroraSimulation initialize! populates workspace" begin
    mktempdir() do savedir
        altitude_lims = [100, 400]
        θ_lims = 180:-45:0
        E_max = 100
        B_angle_to_zenith = 13

        msis_file = find_msis_file(; verbose=false)
        iri_file = find_iri_file(; verbose=false)

        model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)
        flux = InputFlux(FlatSpectrum(1.0; E_min=50.0), SmoothOnset(0.0, 0.05);
                         beams=1:2, z_source=500.0)

        sim = AuroraSimulation(model, flux, savedir;
                               mode=TimeDependentMode(duration = 0.1, dt = 0.01,
                                                      CFL_number = 128, n_loop = 2))

        @test !sim.workspace.initialized
        @test sim.time isa RefinedTimeGrid
        @test sim.time.dt_internal <= sim.time.dt

        initialize!(sim; verbose=false)

        @test sim.workspace.initialized
        @test sim.model.initialized
        n_species = 3
        @test sim.workspace.degradation.secondary_e_flux isa NTuple{n_species, Matrix{Float64}}
        @test sim.workspace.degradation.primary_e_spectrum isa NTuple{n_species, Vector{Float64}}
        @test all(!isempty(sp.cascading_data.E_edges) for sp in sim.model.species)
        @test size(sim.workspace.Ie, 2) == sim.time.n_t_per_loop
        @test size(sim.workspace.Ie_top, 2) == length(sim.time.t)
    end
end

@testitem "Failed rebuild leaves workspace uninitialized" begin
    mktempdir() do savedir
        msis_file = find_msis_file(; verbose=false)
        iri_file = find_iri_file(; verbose=false)

        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
        flux = InputFlux(FlatSpectrum(1e-2; E_min=50.0); beams=1:2, z_source=250.0)
        sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode())
        initialize!(sim; verbose=false)
        @test sim.workspace.initialized

        # Corrupt the mutable beam-index vector to force input-flux construction to fail after
        # initialize! has invalidated the existing workspace.
        empty!(flux.beams)
        @test_throws ErrorException initialize!(sim; verbose=false)
        @test !sim.workspace.initialized
        @test AURORA.needs_initialization(sim)
    end
end

@testitem "AuroraSimulation run! auto-initializes" begin
    mktempdir() do savedir
        altitude_lims = [100, 400]
        θ_lims = 180:-45:0
        E_max = 100
        B_angle_to_zenith = 13
        msis_file = find_msis_file(; verbose=false)
        iri_file = find_iri_file(; verbose=false)

        model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)
        flux = InputFlux(FlatSpectrum(1.0; E_min=50.0); beams=1:2)
        sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode())

        @test !sim.workspace.initialized

        run!(sim; verbose=false)

        @test sim.workspace.initialized
    end
end

@testitem "Multi-step SS: saved t_run matches time grid" begin
    using NCDatasets
    mktempdir() do savedir
        altitude_lims = [100, 200]
        θ_lims = 180:-90:0
        E_max = 100
        B_angle_to_zenith = 13

        msis_file = find_msis_file(; verbose=false)
        iri_file = find_iri_file(; verbose=false)

        model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)
        flux = InputFlux(FlatSpectrum(1.0; E_min=50.0), SinusoidalFlickering(5.0); beams=1:2)
        sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode(duration = 0.04, dt = 0.01))

        run!(sim; verbose=false)

        NCDataset(joinpath(savedir, "simulation_data.nc"), "r") do ds
            t_run = Array(ds["time"])
            expected_t = collect(sim.time.t)

            # t_run must span the full time axis, not be a scalar 1
            @test length(t_run) == length(expected_t)
            @test t_run ≈ expected_t

            # Ie time dimension must also match
            @test size(ds["Ie"], 3) == length(expected_t)
        end
    end
end

@testitem "SteadyStateMode() → SingleStepConfig" begin
    mktempdir() do savedir
        msis_file = find_msis_file(; verbose=false)
        iri_file  = find_iri_file(; verbose=false)
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

@testitem "NeutralSpecies density_source types" begin
    msis_file = find_msis_file(; verbose=false)
    iri_file  = find_iri_file(; verbose=false)
    model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)

    # Default model from an MSIS file: read_msis_file reads the file eagerly and yields a
    # DensityProfile (carrying a provenance origin); density is empty before initialize!
    for sp in model.species
        @test sp.density_source isa DensityProfile
        @test occursin("MSIS file", sp.density_source.origin)
    end
    @test isempty(model.species[1].density)

    # After initialize! density is populated
    initialize!(model; verbose=false)
    ag = model.altitude_grid
    @test !isempty(model.species[1].density)

    # A user DensityProfile built on the file's native grid reproduces the default density
    raw = AURORA.load_msis(msis_file)
    vd  = DensityProfile(raw.data.height_km .* 1e3, raw.data.N2; origin="manual")
    model_vd = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
    model_vd.species[1].density_source = vd
    initialize!(model_vd; verbose=false)
    @test model_vd.species[1].density_source isa DensityProfile
    @test model_vd.species[1].density ≈ model.species[1].density rtol=1e-6

    # A @law source is accepted, and density remains empty until initialize!
    flat_profile = @law h -> fill(1e15, length(h))
    sp_fn = AURORA.N2Species(flat_profile)
    @test sp_fn.density_source isa ExprLaw
    @test sp_fn.density_source === flat_profile
    @test isempty(sp_fn.density)

    # A bare anonymous law is rejected to ensure reproducibility
    @test_throws ArgumentError AURORA.N2Species(h -> fill(1e15, length(h)))
end

@testitem "NeutralAtmosphere as a model atmosphere" begin
    msis_file = find_msis_file(; verbose=false)
    iri_file  = find_iri_file(; verbose=false)

    neutrals = read_msis_file(msis_file)
    @test neutrals isa NeutralAtmosphere
    @test all(haskey(neutrals, s) for s in (:N2, :O2, :O))
    @test neutrals[:N2] isa DensityProfile
    @test_throws ArgumentError neutrals[:XX]

    # A NeutralAtmosphere is accepted wherever an MSIS path is, and gives the same densities
    model_path = AuroraModel((100, 400), 180:-30:0, 100, msis_file, iri_file)
    model_prof = AuroraModel((100, 400), 180:-30:0, 100, neutrals, iri_file)
    initialize!(model_path)
    initialize!(model_prof)
    for name in (:N2, :O2, :O)
        @test model_prof.species[name].density ≈ model_path.species[name].density rtol=1e-12
    end
end

@testitem "AuroraModel species support Symbol indexing" begin
    msis_file = find_msis_file(; verbose=false)
    iri_file  = find_iri_file(; verbose=false)
    model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)

    @test model.species[:N2] === model.species[1]
    @test model.species[:O2] === model.species[2]
    @test model.species[:O] === model.species[3]
    @test_throws KeyError model.species[:NO]
end

@testitem "Species Symbol indexing rejects duplicate names" begin
    msis_file = find_msis_file(; verbose=false)
    iri_file  = find_iri_file(; verbose=false)
    model = AuroraModel([100, 200], 180:-90:0, 100, nothing, iri_file, 0;
                        species = (N2Species(msis_file), N2Species(msis_file)))

    @test_throws ArgumentError model.species[:N2]
end

@testitem "AuroraModel requires exactly one of neutrals or species" begin
    msis_file = find_msis_file(; verbose=false)
    iri_file  = find_iri_file(; verbose=false)

    # Neither given: the model has no way to build the default species.
    @test_throws "either a neutral atmosphere" AuroraModel(
        [100, 200], 180:-90:0, 100, nothing, iri_file, 0)

    # Both given: neutrals would be silently ignored.
    @test_throws "not used when `species`" AuroraModel(
        [100, 200], 180:-90:0, 100, msis_file, iri_file, 0;
        species = (N2Species(msis_file),))

    # neutrals = nothing with an explicit species tuple constructs fine.
    model = AuroraModel([100, 200], 180:-90:0, 100, nothing, iri_file, 0;
                        species = (N2Species(msis_file),))
    @test model isa AuroraModel
end

@testitem "AuroraModel is uninitialized before initialize!" begin
    msis_file = find_msis_file(; verbose=false)
    iri_file  = find_iri_file(; verbose=false)
    model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)

    @test !model.initialized
    @test isempty(model.scattering.θ_scatter)
    @test isempty(model.species[1].density)
end

@testitem "initialize!(model) interception window" begin
    mktempdir() do savedir
        msis_file = find_msis_file(; verbose=false)
        iri_file  = find_iri_file(; verbose=false)
        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)

        flat_n2 = @law h -> fill(1e18, length(h))
        model.species[:N2].density_source = flat_n2

        flux = InputFlux(FlatSpectrum(1.0; E_min=50.0); beams=1:2)
        sim  = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode())

        @test !sim.model.initialized
        run!(sim; verbose=false)
        @test sim.model.initialized

        n2_density = sim.model.species[:N2].density
        @test !isempty(n2_density)
        @test n2_density[1] ≈ 1e18   # bottom of the grid, unaffected by boundary taper
    end
end

@testitem "AuroraModel with 2 species: run! succeeds" begin
    mktempdir() do savedir
        msis_file = find_msis_file(; verbose=false)
        iri_file  = find_iri_file(; verbose=false)

        model = AuroraModel([100, 200], 180:-90:0, 100, nothing, iri_file, 0;
                            species = (O2Species(msis_file), OSpecies(msis_file)))
        flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim  = AuroraSimulation(model, flux, savedir; mode = SteadyStateMode())
        run!(sim; verbose=false)

        @test sim.model.initialized
        @test length(sim.model.species) == 2
        @test sim.workspace.degradation.secondary_e_flux isa NTuple{2, Matrix{Float64}}
    end
end

@testitem "Custom 4th species from a channel table: run! succeeds" begin
    import TOML
    mktempdir() do savedir
        msis_file = find_msis_file(; verbose=false)
        iri_file  = find_iri_file(; verbose=false)

        custom_law = @law (E_s, E_p) -> 1.0 / (11.4^2 + E_s^2)
        # Invented gas: one excitation, one single and one double ionization channel
        channels = [AURORA.CollisionChannel("exc", @law(E -> 1e-21 .* (E .> 8.0)), 8.0, 0),
                    AURORA.CollisionChannel("ion", @law(E -> 5e-22 .* (E .> 20.0)), 20.0, 1;
                                            source = "invented"),
                    AURORA.CollisionChannel("dion", @law(E -> 1e-22 .* (E .> 35.0)), 35.0, 2)]
        elastic = @law(E -> fill(1e-20, length(E)))
        # Reuse phase function from N₂
        custom_sp = AURORA.NeutralSpecies(:CustomGas, @law(h -> fill(1e18, length(h)));
                                          elastic_cross_section = elastic,
                                          channels,
                                          secondary_law         = custom_law,
                                          phase_fcn_generator   = AURORA.phase_fcn_N2)

        model = AuroraModel([100, 200], 180:-90:0, 100, nothing, iri_file, 0;
                            species = (N2Species(msis_file), O2Species(msis_file),
                                       OSpecies(msis_file), custom_sp))

        flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim  = AuroraSimulation(model, flux, savedir; mode = SteadyStateMode())
        run!(sim; verbose=false)

        @test sim.model.initialized
        @test length(sim.model.species) == 4
        @test sim.workspace.degradation.secondary_e_flux isa NTuple{4, Matrix{Float64}}
        @test sim.model.species[end].name == :CustomGas
        @test !isempty(sim.model.species[end].density)
        @test AURORA.channel_names(sim.model.species[end]) == ["exc", "ion", "dion"]
        @test sim.model.species[end].excitation_levels == [0.0 0.0; 8.0 0.0; 20.0 1.0; 35.0 2.0]
        @test sim.model.species[end].cascading_spec.ionization_thresholds == [20.0, 35.0]
        @test sim.model.species[end].cascading_spec.n_secondaries == [1, 2]

        # The channel tables are written to inputs/collision_channels.toml
        channels_toml = TOML.parsefile(joinpath(savedir, "inputs", "collision_channels.toml"))
        @test [c["name"] for c in channels_toml["CustomGas"]["channels"]] ==
              ["exc", "ion", "dion"]
        @test channels_toml["CustomGas"]["channels"][2]["n_secondaries"] == 1
        @test channels_toml["CustomGas"]["channels"][2]["source"] == "invented"
        @test [c["energy_loss_eV"] for c in channels_toml["N2"]["channels"]][end] == 42.0
    end
end

@testitem "Custom phase function via interception window: run! succeeds" begin
    mktempdir() do savedir
        msis_file = find_msis_file(; verbose=false)
        iri_file  = find_iri_file(; verbose=false)

        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)

        custom_generator = @law (θ, E) -> AURORA.phase_fcn_N2(θ, E)
        model.species[1].phase_fcn_generator = custom_generator

        flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim  = AuroraSimulation(model, flux, savedir; mode = SteadyStateMode())
        run!(sim; verbose=false)

        @test sim.model.initialized
        @test model.species[1].phase_fcn_generator === custom_generator
    end
end

@testitem "Altitude grid swap: initialize!(model) rebuilds s_field and ionosphere" begin
    msis_file = find_msis_file(; verbose=false)
    iri_file  = find_iri_file(; verbose=false)

    model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
    initialize!(model; verbose=false)
    old_n_z = model.altitude_grid.n
    @test length(model.s_field)               == old_n_z
    @test length(model.ionosphere.ne)         == old_n_z
    @test length(model.species[1].density)    == old_n_z

    model.altitude_grid = AltitudeGrid(100, 300)
    initialize!(model; verbose=false)

    new_n_z = model.altitude_grid.n
    @test new_n_z > old_n_z
    @test length(model.s_field)               == new_n_z
    @test length(model.ionosphere.ne)         == new_n_z
    @test length(model.species[1].density)    == new_n_z
end

@testitem "Altitude grid swap: run! after initialize!(model) succeeds" begin
    mktempdir() do savedir
        msis_file = find_msis_file(; verbose=false)
        iri_file  = find_iri_file(; verbose=false)

        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
        flux  = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        output = AuroraOutputManager(savedir; overwrite=true)
        sim   = AuroraSimulation(model, flux, output; mode = SteadyStateMode())
        run!(sim; verbose=false)

        model.altitude_grid = AltitudeGrid(100, 300)
        initialize!(model; verbose=false)   # recomputes s_field, ionosphere, species
        initialize!(sim; verbose=false)     # rebuilds workspace for new grid dimensions
        run!(sim)

        @test sim.model.initialized
        @test length(sim.model.s_field) == model.altitude_grid.n
    end
end

@testitem "Reassigning model inputs invalidates the model" begin
    msis_file = find_msis_file(; verbose=false)
    iri_file  = find_iri_file(; verbose=false)

    model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
    initialize!(model; verbose=false)
    @test model.initialized

    # Each geometry reassignment must flip `initialized` back to false.
    model.altitude_grid = AltitudeGrid(100, 300)
    @test !model.initialized

    initialize!(model; verbose=false)
    model.energy_grid = EnergyGrid(200)
    @test !model.initialized

    initialize!(model; verbose=false)
    model.B_angle_to_zenith = 20
    @test !model.initialized

    initialize!(model; verbose=false)
    replacement_species = deepcopy(model.species)
    model.species = replacement_species
    @test !model.initialized
end

@testitem "Grid change then bare run! auto-reinitializes (no manual init)" begin
    mktempdir() do savedir
        msis_file = find_msis_file(; verbose=false)
        iri_file  = find_iri_file(; verbose=false)

        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
        flux  = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        output = AuroraOutputManager(savedir; overwrite=true)
        sim   = AuroraSimulation(model, flux, output; mode = SteadyStateMode())
        run!(sim; verbose=false)

        # Change the grid and call run! directly — no initialize!(model)/initialize!(sim).
        model.altitude_grid = AltitudeGrid(100, 300)
        run!(sim; verbose=false)

        @test sim.model.initialized
        @test length(sim.model.s_field) == model.altitude_grid.n
        @test size(sim.workspace.Ie, 1) ÷ length(model.pitch_angle_grid.μ_center) == model.altitude_grid.n
    end
end

@testitem "AuroraOutputManager compress kwarg" begin
    # true/false/integer conversion and out-of-range guard
    @test AuroraOutputManager("x"; compress=true).deflatelevel  == 4
    @test AuroraOutputManager("x"; compress=false).deflatelevel == 0
    @test AuroraOutputManager("x"; compress=6).deflatelevel     == 6
    @test AuroraOutputManager("x"; compress=0).deflatelevel     == 0
    @test_throws ArgumentError AuroraOutputManager("x"; compress=10)
    @test_throws ArgumentError AuroraOutputManager("x"; compress=-1)
end

@testitem "AuroraOutputManager savedir normalization" begin
    for savedir in ("", " ", "\t", "\n", " \t\n")
        fallback = AuroraOutputManager(savedir).savedir
        @test dirname(fallback) == "backup"
        @test occursin(r"^\d{8}-\d{4}$", basename(fallback))
    end

    @test AuroraOutputManager("path with spaces").savedir == "path with spaces"
end

@testitem "Higher compress level produces smaller simulation_data.nc" begin
    msis_file = find_msis_file(; verbose=false)
    iri_file  = find_iri_file(; verbose=false)
    model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 13)
    flux  = InputFlux(FlatSpectrum(1.0; E_min=50.0); beams=1:2)

    # Use a multi-step run so the Ie array is big enough for the compression to have an effect
    mode = SteadyStateMode(duration=0.5, dt=0.01)

    size_lo = mktempdir() do dir
        sim = AuroraSimulation(model, flux, AuroraOutputManager(dir; compress=false); mode)
        run!(sim; verbose=false)
        filesize(joinpath(dir, "simulation_data.nc"))
    end

    size_hi = mktempdir() do dir
        sim = AuroraSimulation(model, flux, AuroraOutputManager(dir; compress=9); mode)
        run!(sim; verbose=false)
        filesize(joinpath(dir, "simulation_data.nc"))
    end

    @test size_hi < size_lo
end

@testitem "Energy grid change rebuilds sim.time and workspace (TimeDependent)" begin
    mktempdir() do savedir
        msis_file = find_msis_file(; verbose=false)
        iri_file  = find_iri_file(; verbose=false)

        model = AuroraModel([100, 200], 180:-90:0, 80, msis_file, iri_file, 0)
        flux  = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        output = AuroraOutputManager(savedir; overwrite=true)
        sim   = AuroraSimulation(model, flux, output;
                                 mode = TimeDependentMode(duration=0.02, dt=0.01,
                                                          CFL_number=128, n_loop=1))
        run!(sim; verbose=false)
        @test size(sim.workspace.Ie, 3) == model.energy_grid.n

        # Larger energy grid → more energy bins AND a different CFL-refined time grid.
        model.energy_grid = EnergyGrid(200)
        run!(sim; verbose=false)

        @test sim.model.initialized
        @test size(sim.workspace.Ie, 3) == model.energy_grid.n
        @test sim.time isa AURORA.RefinedTimeGrid
    end
end

@testitem "Law enforcement: bare lambdas and captured locals rejected" begin
    msis_file = find_msis_file(; verbose=false)
    iri_file  = find_iri_file(; verbose=false)

    default_kwargs = (; elastic_cross_section = AURORA.default_elastic_cross_section(:N2),
                        channels              = AURORA.default_channels(:N2),
                        secondary_law         = AURORA.default_secondary_law(:N2),
                        phase_fcn_generator   = AURORA.phase_fcn_N2)

    # Bare anonymous functions are rejected
    @test_throws ArgumentError AURORA.CascadingSpec("X", [1.0], (a, b) -> a)
    @test_throws ArgumentError AURORA.N2Species(h -> fill(1e15, length(h)))
    @test_throws ArgumentError AURORA.NeutralSpecies(:G, read_msis_file(msis_file)[:N2];
                                   default_kwargs..., phase_fcn_generator = (θ, E) -> θ)
    @test_throws ArgumentError AURORA.NeutralSpecies(:G, read_msis_file(msis_file)[:N2];
                                   default_kwargs..., elastic_cross_section = E -> E)

    # A @law that closes over a local variable is rejected (its source can't be rebuilt)
    @test_throws ArgumentError (let n0 = 1e18
        @law h -> fill(n0, length(h))
    end)

    # A non-callable object is rejected too — the realistic mistake of assigning a whole
    # NeutralAtmosphere as density_source instead of indexing it (neutrals[:N2])
    @test_throws "must be callable" AURORA.NeutralSpecies(:G, read_msis_file(msis_file);
                                   default_kwargs...)

    # @law, functors and named functions are all accepted
    @test (@law h -> fill(1e15, length(h))) isa ExprLaw
    sp = AURORA.N2Species(read_msis_file(msis_file)[:N2])
    @test sp.density_source isa DensityProfile        # eager file read → DensityProfile
    @test sp.phase_fcn_generator === phase_fcn_N2    # named function
end

@testitem "The positional NeutralSpecies constructor enforces reproducibility" begin
    empty_mat = Matrix{Float64}(undef, 0, 0)
    spec      = AURORA.default_cascading_spec(:N2)
    cache     = AURORA.SpeciesCascadingCache(spec)
    channels  = AURORA.default_channels(:N2)
    elastic   = AURORA.default_elastic_cross_section(:N2)
    law       = AURORA.default_secondary_law(:N2)

    positional(density_source, elastic_cross_section, secondary_law, phase_fcn_generator) =
        AURORA.NeutralSpecies(:G, density_source, Float64[], elastic_cross_section, channels,
                              secondary_law, phase_fcn_generator, (empty_mat, copy(empty_mat)),
                              copy(empty_mat), copy(empty_mat), spec, cache)

    @test_throws "bare anonymous function" positional(h -> fill(1e18, length(h)), elastic,
                                                      law, AURORA.phase_fcn_N2)
    @test_throws "bare anonymous function" positional(@law(h -> fill(1e18, length(h))),
                                                      E -> E, law, AURORA.phase_fcn_N2)
    @test_throws "bare anonymous function" positional(@law(h -> fill(1e18, length(h))),
                                                      elastic, (E_s, E_p) -> E_s,
                                                      AURORA.phase_fcn_N2)
    @test_throws "bare anonymous function" positional(@law(h -> fill(1e18, length(h))),
                                                      elastic, law, (θ, E) -> θ)

    sp = positional(@law(h -> fill(1e18, length(h))), elastic, law, AURORA.phase_fcn_N2)
    @test sp.name === :G
end

@testitem "@law density round-trips through physics_state.jld2" begin
    using JLD2
    mktempdir() do savedir
        msis_file = find_msis_file(; verbose=false)
        iri_file  = find_iri_file(; verbose=false)

        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
        model.species[:N2].density_source = @law h -> fill(1e18, length(h))
        flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim  = AuroraSimulation(model, flux, savedir; mode = SteadyStateMode())
        run!(sim; verbose=false)

        savefile = joinpath(savedir, "inputs", "physics_state.jld2")
        model2 = JLD2.load(savefile, "model")

        prof = model2.species[:N2].density_source
        @test prof isa ExprLaw
        @test prof.src == model.species[:N2].density_source.src
        # Reconstructed law is callable in this same scope (relies on invokelatest)
        h = model2.altitude_grid.h
        @test prof(h) == fill(1e18, length(h))
        # Reloaded model re-initializes from the reconstructed law
        initialize!(model2; verbose=false)
        @test model2.species[:N2].density[1] ≈ 1e18
    end
end

@testitem "ElectronProfile round-trips through physics_state.jld2" begin
    using JLD2
    mktempdir() do savedir
        msis_file = find_msis_file(; verbose=false)
        iri_file  = find_iri_file(; verbose=false)

        # Build the electron background as an ElectronProfile (no file path stored on the model)
        electron = read_iri_file(iri_file)
        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, electron, 0)
        @test model.ionosphere.electron_source isa ElectronProfile
        flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim  = AuroraSimulation(model, flux, savedir; mode = SteadyStateMode())
        run!(sim; verbose=false)

        model2 = JLD2.load(joinpath(savedir, "inputs", "physics_state.jld2"), "model")
        es = model2.ionosphere.electron_source
        @test es isa ElectronProfile
        @test es.origin == electron.origin
        # Reloaded model re-samples electrons from the stored profile (no external file)
        initialize!(model2; verbose=false)
        @test model2.ionosphere.ne ≈ model.ionosphere.ne rtol=1e-9
        @test model2.ionosphere.Te ≈ model.ionosphere.Te rtol=1e-9
    end
end

@testitem "Default model laws round-trip through physics_state.jld2" begin
    using JLD2, NCDatasets
    mktempdir() do savedir
        msis_file = find_msis_file(; verbose=false)
        iri_file  = find_iri_file(; verbose=false)

        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
        flux  = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim   = AuroraSimulation(model, flux, savedir; mode = SteadyStateMode())
        run!(sim; verbose=false)
        Ie1 = NCDataset(joinpath(savedir, "simulation_data.nc"), "r") do ds
            Array(ds["Ie"])
        end

        model2 = JLD2.load(joinpath(savedir, "inputs", "physics_state.jld2"), "model")

        # Prove laws round-tripped by running a full simulation from the reloaded model
        # and checking that results are identical
        mktempdir() do savedir2
            sim2 = AuroraSimulation(model2, flux, savedir2; mode = SteadyStateMode())
            run!(sim2; verbose=false)
            Ie2 = NCDataset(joinpath(savedir2, "simulation_data.nc"), "r") do ds
                Array(ds["Ie"])
            end
            @test all(Ie2 .≈ Ie1)
        end
    end
end

@testitem "Channel tables round-trip through physics_state.jld2" begin
    using JLD2
    mktempdir() do savedir
        msis_file = find_msis_file(; verbose=false)
        iri_file  = find_iri_file(; verbose=false)

        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
        sp = model.species[:O2]
        # A channel whose cross section is an @law, and a renamed built-in channel
        push!(sp.channels, AURORA.CollisionChannel("mystate",
                                                   @law(E -> fill(1e-21, length(E))),
                                                   9.0, 0; source = "invented"))
        sp.channels[1] = AURORA.CollisionChannel(sp.channels[1]; name = "renamed")

        flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim  = AuroraSimulation(model, flux, savedir; mode = SteadyStateMode())
        run!(sim; verbose=false)

        model2 = JLD2.load(joinpath(savedir, "inputs", "physics_state.jld2"), "model")
        sp2 = model2.species[:O2]
        @test AURORA.channel_names(sp2) == AURORA.channel_names(sp)
        @test AURORA.channel(sp2, "mystate").cross_section isa ExprLaw
        @test AURORA.channel(sp2, "mystate").source == "invented"

        # Wipe the saved matrices so the comparison can only pass if initialize! rebuilds
        # them from the reloaded channel table
        setfield!(sp2, :cross_sections, Matrix{Float64}(undef, 0, 0))
        setfield!(sp2, :excitation_levels, Matrix{Float64}(undef, 0, 0))
        @test_throws "derived from its `elastic_cross_section`" sp2.cross_sections = sp.cross_sections
        initialize!(model2; verbose=false)
        @test sp2.cross_sections == sp.cross_sections
        @test sp2.excitation_levels == sp.excitation_levels

        # ... and that the rebuild tracks a further edit of the reloaded table
        i = findfirst(c -> c.name == "mystate", sp2.channels)
        sp2.channels[i] = AURORA.CollisionChannel(sp2.channels[i]; energy_loss = 3.0)
        initialize!(model2; verbose=false)
        @test sp2.excitation_levels[i + 1, 1] == 3.0
        @test sp2.excitation_levels != sp.excitation_levels
    end
end

@testitem "Editing a channel before initialize! changes the derived data" begin
    msis_file = find_msis_file(; verbose=false)
    iri_file  = find_iri_file(; verbose=false)

    model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
    sp = model.species[:N2]
    i = findfirst(c -> c.name == "a3sup", sp.channels)
    sp.channels[i] = AURORA.CollisionChannel(sp.channels[i]; energy_loss = 14.0)
    filter!(c -> c.name != "ddion", sp.channels)

    initialize!(model; verbose=false)

    @test sp.excitation_levels[i + 1, 1] == 14.0
    @test size(sp.excitation_levels, 1) == length(sp.channels) + 1
    @test size(sp.cross_sections, 1) == length(sp.channels) + 1
    @test 42.0 ∉ sp.cascading_spec.ionization_thresholds
    @test sp.cascading_spec.n_secondaries == [1, 1, 1, 1]
end

@testitem "Editing a channel after initialize! takes effect on the next initialize!" begin
    mktempdir() do savedir
        msis_file = find_msis_file(; verbose=false)
        iri_file  = find_iri_file(; verbose=false)

        model = AuroraModel([100, 200], 180:-90:0, 100, msis_file, iri_file, 0)
        initialize!(model; verbose=false)
        sp = model.species[:N2]
        spec_before  = sp.cascading_spec
        cache_before = sp.cascading_data

        # An excitation channel: the cascading spec is untouched, so the cache object stays
        i = findfirst(c -> c.name == "a3sup", sp.channels)
        sp.channels[i] = AURORA.CollisionChannel(sp.channels[i]; energy_loss = 6.0)
        initialize!(model; verbose=false)
        @test sp.excitation_levels[i + 1, 1] == 6.0
        @test sp.cross_sections[i + 1, :] == AURORA.e_N2a3sup(model.energy_grid.E_centers)
        @test sp.cascading_data === cache_before
        @test sp.cascading_spec.ionization_thresholds == spec_before.ionization_thresholds

        # An ionizing channel: the spec changes, so the cache is replaced
        filter!(c -> c.name != "dion", sp.channels)
        initialize!(model; verbose=false)
        @test sp.cascading_spec.ionization_thresholds == [15.581, 16.73, 18.75, 42.0]
        @test sp.cascading_spec.n_secondaries == [1, 1, 1, 2]
        @test sp.cascading_data !== cache_before
        @test size(sp.excitation_levels, 1) == length(sp.channels) + 1

        # The edited model still runs
        flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim  = AuroraSimulation(model, flux, savedir; mode = SteadyStateMode())
        run!(sim; verbose=false)
        @test sim.model.initialized
    end
end
