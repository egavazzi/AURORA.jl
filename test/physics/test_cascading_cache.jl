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

@testitem "Custom CascadingSpec produces valid transfer matrices" begin
    flat_law = (E_s, E_p) -> 1.0
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
    @test all(cache.secondary_transfer_matrix[1:i_threshold, :, :] .== 0)
end

@testitem "Custom NeutralSpecies with custom CascadingSpec: spectra are accessible" begin
    msis_file = find_msis_file(; verbose=false)

    custom_law  = (E_s, E_p) -> 1.0 / (11.4^2 + E_s^2)
    custom_spec = AURORA.CascadingSpec("N2variant", [15.581, 16.73, 18.75], custom_law)

    sp = AURORA.NeutralSpecies(:N2, AURORA.MSISDensity(msis_file, :N2);
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
