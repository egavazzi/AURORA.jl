@testitem "SpeciesCascadingCache provides spectra accessors" begin
    energy_grid = AURORA.EnergyGrid(100)
    cache = AURORA.SpeciesCascadingCache(AURORA.DefaultCascadingSpecN2())

    AURORA.load_or_compute_cascading!(cache, energy_grid;
                                      policy=AURORA.CachePolicy(force_recompute=true, save_cache=false),
                                      verbose=false)
    i_primary = searchsortedlast(cache.E_edges, 40.0)
    secondary = AURORA.secondary_spectrum(cache, i_primary, 15.581)
    primary = AURORA.primary_spectrum(cache, i_primary, 15.581)

    @test length(secondary) == energy_grid.n
    @test length(primary) == energy_grid.n
    @test sum(secondary) > 0
    @test sum(primary) >= 0
    @test !isempty(cache.E_edges)
    @test !isempty(cache.ionization_thresholds)
    @test size(cache.secondary_transfer_matrix, 1) == energy_grid.n
    @test size(cache.secondary_transfer_matrix, 2) == energy_grid.n
end

@testitem "Cascading cache lifecycle" begin
    using JLD2: jldopen

    cache_files(dir) = isdir(dir) ? filter(name -> endswith(name, ".jld2"), readdir(dir)) : String[]
    compatible_cache_count(dir) = count(cache_files(dir)) do filename
        filepath = joinpath(dir, filename)
        jldopen(filepath, "r") do file
            string(file["version_AURORA"]) == AURORA.cache_version_string()
        end
    end

    cache_root = mktempdir()
    energy_grid = AURORA.EnergyGrid(60)

    save_policy = AURORA.CachePolicy(force_recompute=true, save_cache=true, cache_root=cache_root)
    load_policy = AURORA.CachePolicy(cache_root=cache_root)
    skip_save_policy = AURORA.CachePolicy(force_recompute=true, save_cache=false,
                                          cache_root=joinpath(cache_root, "skip_save"))

    n2_cache = AURORA.SpeciesCascadingCache(AURORA.DefaultCascadingSpecN2())
    o2_cache = AURORA.SpeciesCascadingCache(AURORA.DefaultCascadingSpecO2())
    o_cache  = AURORA.SpeciesCascadingCache(AURORA.DefaultCascadingSpecO())
    AURORA.load_or_compute_cascading!(n2_cache, energy_grid; policy=save_policy, verbose=false)
    AURORA.load_or_compute_cascading!(o2_cache, energy_grid; policy=save_policy, verbose=false)
    AURORA.load_or_compute_cascading!(o_cache,  energy_grid; policy=save_policy, verbose=false)

    n2_dir = joinpath(cache_root, "e_cascading", "N2")
    o2_dir = joinpath(cache_root, "e_cascading", "O2")
    o_dir  = joinpath(cache_root, "e_cascading", "O")
    @test length(cache_files(n2_dir)) == 1
    @test length(cache_files(o2_dir)) == 1
    @test length(cache_files(o_dir))  == 1

    loaded_cache = AURORA.SpeciesCascadingCache(AURORA.DefaultCascadingSpecN2())
    AURORA.load_or_compute_cascading!(loaded_cache, energy_grid; policy=load_policy, verbose=false)
    @test loaded_cache.E_edges == n2_cache.E_edges
    @test loaded_cache.ionization_thresholds == n2_cache.ionization_thresholds

    n2_file = joinpath(n2_dir, only(cache_files(n2_dir)))
    payload = jldopen(n2_file, "r") do file
        (
            Q_primary   = file["Q_primary"],
            Q_secondary = file["Q_secondary"],
            E_edges     = file["E_edges"],
            E_ionizations = file["E_ionizations"],
        )
    end
    rm(n2_file; force=true)
    jldopen(n2_file, "w") do file
        file["version_AURORA"] = "0.0.0"
        file["Q_primary"]      = payload.Q_primary
        file["Q_secondary"]    = payload.Q_secondary
        file["E_edges"]        = payload.E_edges
        file["E_ionizations"]  = payload.E_ionizations
    end

    stale_cache = AURORA.SpeciesCascadingCache(AURORA.DefaultCascadingSpecN2())
    AURORA.load_or_compute_cascading!(stale_cache, energy_grid; policy=load_policy, verbose=false)
    @test compatible_cache_count(n2_dir) >= 1

    AURORA.load_or_compute_cascading!(AURORA.SpeciesCascadingCache(AURORA.DefaultCascadingSpecN2()),
                                      energy_grid; policy=skip_save_policy, verbose=false)
    skip_n2_dir = joinpath(cache_root, "skip_save", "e_cascading", "N2")
    @test isempty(cache_files(skip_n2_dir))

    AURORA.clear_cascading_cache!(cache_root=cache_root)
    @test isempty(cache_files(n2_dir))
    @test isempty(cache_files(o2_dir))
    @test isempty(cache_files(o_dir))
end

@testitem "Cascading cache built for another spec is rejected" begin
    cache_files(dir) = isdir(dir) ? filter(name -> endswith(name, ".jld2"), readdir(dir)) : String[]

    cache_root = mktempdir()
    energy_grid = AURORA.EnergyGrid(60)
    n2_dir = joinpath(cache_root, "e_cascading", "N2")

    # A spec sharing the species name (and hence the cache directory) with the default N₂ one,
    # but describing different ionization physics.
    other_law = @law (E_s, E_p) -> 1.0 / (9.0^2 + E_s^2)
    other_spec = AURORA.CascadingSpec("N2", [99.0], other_law)
    other_cache = AURORA.SpeciesCascadingCache(other_spec)
    AURORA.load_or_compute_cascading!(other_cache, energy_grid; verbose = false,
        policy = AURORA.CachePolicy(; force_recompute = true, save_cache = true, cache_root))
    @test length(cache_files(n2_dir)) == 1

    n2_cache = AURORA.SpeciesCascadingCache(AURORA.DefaultCascadingSpecN2())
    AURORA.load_or_compute_cascading!(n2_cache, energy_grid; verbose = false,
                                      policy = AURORA.CachePolicy(; cache_root))

    @test size(n2_cache.primary_transfer_matrix, 3) == 5
    @test size(n2_cache.secondary_transfer_matrix, 3) == 5
    @test n2_cache.ionization_thresholds == AURORA.DefaultCascadingSpecN2().ionization_thresholds
    # A compatible file was saved, so a second request loads from disk.
    found, _ = AURORA.find_cascading_cache(AURORA.DefaultCascadingSpecN2(), energy_grid.E_edges;
                                           verbose = false,
                                           policy = AURORA.CachePolicy(; cache_root))
    @test found

    # Same thresholds and secondary counts, different law: also incompatible.
    law_only_spec = AURORA.CascadingSpec("N2", AURORA.DefaultCascadingSpecN2().ionization_thresholds,
                                         other_law;
                                         n_secondaries = AURORA.DefaultCascadingSpecN2().n_secondaries)
    found, _ = AURORA.find_cascading_cache(law_only_spec, energy_grid.E_edges; verbose = false,
                                           policy = AURORA.CachePolicy(; cache_root))
    @test !found

    # Same thresholds and law, different secondary counts: also incompatible.
    n_sec_spec = AURORA.CascadingSpec("N2", AURORA.DefaultCascadingSpecN2().ionization_thresholds,
                                      AURORA.DefaultCascadingSpecN2().secondary_law;
                                      n_secondaries = [1, 1, 1, 2, 2])
    found, _ = AURORA.find_cascading_cache(n_sec_spec, energy_grid.E_edges; verbose = false,
                                           policy = AURORA.CachePolicy(; cache_root))
    @test !found
end

@testitem "Cascading cache tracks the secondary law's parameters" begin
    struct ScaledLorentzian
        width::Float64
    end
    (law::ScaledLorentzian)(E_s, E_p) = 1.0 / (law.width^2 + E_s^2)

    cache_root = mktempdir()
    energy_grid = AURORA.EnergyGrid(60)

    spec_a = AURORA.CascadingSpec("Functor", [15.581], ScaledLorentzian(11.4))
    spec_b = AURORA.CascadingSpec("Functor", [15.581], ScaledLorentzian(15.2))
    @test AURORA.is_fingerprintable(spec_a.secondary_law)
    @test AURORA.law_fingerprint(spec_a.secondary_law) !=
          AURORA.law_fingerprint(spec_b.secondary_law)

    cache_a = AURORA.SpeciesCascadingCache(spec_a)
    AURORA.load_or_compute_cascading!(cache_a, energy_grid; verbose = false,
        policy = AURORA.CachePolicy(; force_recompute = true, save_cache = true, cache_root))

    load_policy = AURORA.CachePolicy(; cache_root)
    found_a, _ = AURORA.find_cascading_cache(spec_a, energy_grid.E_edges; verbose = false,
                                             policy = load_policy)
    found_b, _ = AURORA.find_cascading_cache(spec_b, energy_grid.E_edges; verbose = false,
                                             policy = load_policy)
    @test found_a
    @test !found_b

    # A plain named function carries neither a source nor parameters, so it cannot be told
    # apart from another definition of the same name and is never cached.
    flat_law(E_s, E_p) = 1.0
    @test !AURORA.is_fingerprintable(flat_law)
    @test_throws ArgumentError AURORA.law_fingerprint(flat_law)

    spec_fn = AURORA.CascadingSpec("NamedFunction", [15.581], flat_law)
    cache_fn = AURORA.SpeciesCascadingCache(spec_fn)
    AURORA.load_or_compute_cascading!(cache_fn, energy_grid; verbose = false,
        policy = AURORA.CachePolicy(; cache_root))
    @test size(cache_fn.primary_transfer_matrix, 3) == 1
    @test !isdir(joinpath(cache_root, "e_cascading", "NamedFunction"))
end

@testitem "Cascading cache missing a required entry is rejected" begin
    using JLD2: jldopen

    cache_root = mktempdir()
    energy_grid = AURORA.EnergyGrid(60)
    spec = AURORA.DefaultCascadingSpecN2()
    cache = AURORA.SpeciesCascadingCache(spec)
    AURORA.load_or_compute_cascading!(cache, energy_grid; verbose = false,
        policy = AURORA.CachePolicy(; force_recompute = true, save_cache = true, cache_root))

    n2_dir = joinpath(cache_root, "e_cascading", "N2")
    n2_file = joinpath(n2_dir, only(filter(f -> endswith(f, ".jld2"), readdir(n2_dir))))
    payload = jldopen(n2_file, "r") do file
        Dict(key => file[key] for key in AURORA.CASCADING_CACHE_KEYS)
    end
    load_policy = AURORA.CachePolicy(; cache_root)

    # Every entry is required, the version string included, so each one alone disqualifies
    # the file. The others are written unchanged, so nothing else can explain the rejection.
    for dropped in AURORA.CASCADING_CACHE_KEYS
        rm(n2_file; force = true)
        jldopen(n2_file, "w") do file
            for (key, value) in payload
                key == dropped || (file[key] = value)
            end
        end
        found, _ = AURORA.find_cascading_cache(spec, energy_grid.E_edges; verbose = false,
                                               policy = load_policy)
        @test !found
    end

    rm(n2_file; force = true)
    jldopen(n2_file, "w") do file
        for (key, value) in payload
            file[key] = value
        end
    end
    found, _ = AURORA.find_cascading_cache(spec, energy_grid.E_edges; verbose = false,
                                           policy = load_policy)
    @test found
end

@testitem "Cascading spectra accessors require an exact threshold" begin
    energy_grid = AURORA.EnergyGrid(100)
    cache = AURORA.SpeciesCascadingCache(AURORA.DefaultCascadingSpecN2())
    AURORA.load_or_compute_cascading!(cache, energy_grid;
                                      policy = AURORA.CachePolicy(force_recompute = true,
                                                                  save_cache = false),
                                      verbose = false)

    @test_throws ArgumentError AURORA.primary_spectrum(cache, 60, 15.6)
    @test_throws ArgumentError AURORA.secondary_spectrum(cache, 60, 15.6)
    @test_throws "available thresholds" AURORA.primary_spectrum(cache, 60, 20.0)

    # The energy-argument methods select the bin that contains the given energy.
    k = searchsortedlast(cache.E_edges, 40.0)
    E_center = (cache.E_edges[k] + cache.E_edges[k + 1]) / 2
    @test AURORA.primary_spectrum(cache, E_center, 15.581) == AURORA.primary_spectrum(cache, k, 15.581)
    @test AURORA.secondary_spectrum(cache, E_center, 15.581) == AURORA.secondary_spectrum(cache, k, 15.581)
end

@testitem "Ionizing levels must agree with the cascading spec" begin
    msis_file = find_msis_file(; verbose = false)
    neutrals = read_msis_file(msis_file)

    sp = AURORA.N2Species(neutrals)
    sp.excitation_levels = AURORA.load_excitation_threshold_for("N2")
    @test AURORA.validate_ionization_channels(sp) === nothing

    i_ion = findfirst(i -> sp.excitation_levels[i, 2] > 0, axes(sp.excitation_levels, 1))

    shifted = copy(sp.excitation_levels)
    shifted[i_ion, 1] += 0.5
    sp.excitation_levels = shifted
    @test_throws ArgumentError AURORA.validate_ionization_channels(sp)

    miscounted = AURORA.load_excitation_threshold_for("N2")
    miscounted[i_ion, 2] = 2
    sp.excitation_levels = miscounted
    @test_throws ArgumentError AURORA.validate_ionization_channels(sp)

    fractional = AURORA.load_excitation_threshold_for("N2")
    fractional[i_ion, 2] = 0.5
    sp.excitation_levels = fractional
    @test_throws ArgumentError AURORA.validate_ionization_channels(sp)
end

@testitem "Custom CascadingSpec produces valid transfer matrices" begin
    flat_law = @law (E_s, E_p) -> 1.0
    custom_spec = AURORA.CascadingSpec("FlatTest", [15.0, 25.0], flat_law)
    cache = AURORA.SpeciesCascadingCache(custom_spec)
    energy_grid = AURORA.EnergyGrid(100)

    AURORA.load_or_compute_cascading!(cache, energy_grid;
                                      policy=AURORA.CachePolicy(save_cache=false),
                                      verbose=false)

    n_E = energy_grid.n
    @test !isempty(cache.E_edges)
    @test length(cache.ionization_thresholds) == 2
    @test size(cache.secondary_transfer_matrix, 1) == n_E
    @test size(cache.secondary_transfer_matrix, 2) == n_E
    @test size(cache.primary_transfer_matrix, 1)   == n_E
    @test size(cache.primary_transfer_matrix, 2)   == n_E
    @test all(cache.secondary_transfer_matrix .>= 0)
    @test all(cache.primary_transfer_matrix   .>= 0)
    # no secondary production below the lowest ionization threshold
    i_threshold = searchsortedlast(cache.E_edges, 15.0)
    @test all(cache.secondary_transfer_matrix[1:(i_threshold - 1), :, :] .== 0)
end

@testitem "Custom NeutralSpecies with custom CascadingSpec: spectra are accessible" begin
    msis_file = find_msis_file(; verbose=false)

    custom_law  = @law (E_s, E_p) -> 1.0 / (11.4^2 + E_s^2)
    custom_spec = AURORA.CascadingSpec("N2variant", [15.581, 16.73, 18.75], custom_law)

    sp = AURORA.NeutralSpecies(:N2, read_msis_file(msis_file)[:N2];
                               cascading_spec      = custom_spec,
                               phase_fcn_generator = AURORA.phase_fcn_N2)

    @test sp.cascading_spec.name == "N2variant"
    @test isempty(sp.density)

    energy_grid = AURORA.EnergyGrid(100)
    AURORA.load_or_compute_cascading!(sp.cascading_data, energy_grid;
                                      policy=AURORA.CachePolicy(save_cache=false),
                                      verbose=false)

    n_E = energy_grid.n
    i_primary = searchsortedlast(sp.cascading_data.E_edges, 40.0)
    secondary = AURORA.secondary_spectrum(sp.cascading_data, i_primary, 15.581)
    primary   = AURORA.primary_spectrum(sp.cascading_data,   i_primary, 15.581)

    @test length(secondary) == n_E
    @test length(primary)   == n_E
    @test sum(secondary) >= 0
    @test all(secondary .>= 0)
    @test all(primary   .>= 0)
end

@testitem "Two quick cascading-cache saves for one species do not collide" begin
    cache_files(dir) = isdir(dir) ? filter(name -> endswith(name, ".jld2"), readdir(dir)) : String[]

    cache_root = mktempdir()
    n2_dir = joinpath(cache_root, "e_cascading", "N2")
    save_policy = AURORA.CachePolicy(; force_recompute = true, save_cache = true, cache_root)

    # Two different grids for the same species, saved back to back: fast enough that the
    # "yyyymmdd-HHMMSS" timestamp in the filename can be identical for both.
    grid_a = AURORA.EnergyGrid(60)
    grid_b = AURORA.EnergyGrid(90)

    cache_a = AURORA.SpeciesCascadingCache(AURORA.DefaultCascadingSpecN2())
    cache_b = AURORA.SpeciesCascadingCache(AURORA.DefaultCascadingSpecN2())
    AURORA.load_or_compute_cascading!(cache_a, grid_a; policy = save_policy, verbose = false)
    AURORA.load_or_compute_cascading!(cache_b, grid_b; policy = save_policy, verbose = false)

    files = cache_files(n2_dir)
    @test length(files) == 2
    @test allunique(files)

    # Both files are intact and independently loadable.
    reloaded_a = AURORA.SpeciesCascadingCache(AURORA.DefaultCascadingSpecN2())
    reloaded_b = AURORA.SpeciesCascadingCache(AURORA.DefaultCascadingSpecN2())
    load_policy_a = AURORA.CachePolicy(; cache_root)
    AURORA.load_or_compute_cascading!(reloaded_a, grid_a; policy = load_policy_a, verbose = false)
    AURORA.load_or_compute_cascading!(reloaded_b, grid_b; policy = load_policy_a, verbose = false)
    @test reloaded_a.E_edges == cache_a.E_edges
    @test reloaded_b.E_edges == cache_b.E_edges
end
