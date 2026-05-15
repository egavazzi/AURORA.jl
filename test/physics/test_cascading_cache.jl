@testitem "CascadingCache provides spectra accessors" begin
    energy_grid = AURORA.EnergyGrid(100)
    cache = AURORA.CascadingCache()

    AURORA.load_or_compute_cascading!(cache, energy_grid)
    i_primary = searchsortedlast(cache[1].E_edges, 40.0)
    secondary = AURORA.secondary_spectrum(cache[1], i_primary, 15.581)
    primary = AURORA.primary_spectrum(cache[1], i_primary, 15.581)

    @test length(cache) == 3
    @test length(secondary) == energy_grid.n
    @test length(primary) == energy_grid.n
    @test sum(secondary) > 0
    @test sum(primary) >= 0
    @test !isempty(cache[1].E_edges)
    @test !isempty(cache[1].ionization_thresholds)
    @test size(cache[1].secondary_transfer_matrix, 1) == energy_grid.n
    @test size(cache[1].secondary_transfer_matrix, 2) == energy_grid.n
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

    cache = AURORA.CascadingCache()
    AURORA.load_or_compute_cascading!(cache, energy_grid; policy=save_policy)

    n2_dir = joinpath(cache_root, "e_cascading", "N2")
    o2_dir = joinpath(cache_root, "e_cascading", "O2")
    o_dir = joinpath(cache_root, "e_cascading", "O")
    @test length(cache_files(n2_dir)) == 1
    @test length(cache_files(o2_dir)) == 1
    @test length(cache_files(o_dir)) == 1

    loaded_cache = AURORA.CascadingCache()
    AURORA.load_or_compute_cascading!(loaded_cache, energy_grid; policy=load_policy)
    @test loaded_cache[1].E_edges == cache[1].E_edges
    @test loaded_cache[1].ionization_thresholds == cache[1].ionization_thresholds

    n2_file = joinpath(n2_dir, only(cache_files(n2_dir)))
    payload = jldopen(n2_file, "r") do file
        (
            Q_primary = file["Q_primary"],
            Q_secondary = file["Q_secondary"],
            E_edges = file["E_edges"],
            E_ionizations = file["E_ionizations"],
        )
    end
    rm(n2_file; force=true)
    jldopen(n2_file, "w") do file
        file["version_AURORA"] = "0.0.0"
        file["Q_primary"] = payload.Q_primary
        file["Q_secondary"] = payload.Q_secondary
        file["E_edges"] = payload.E_edges
        file["E_ionizations"] = payload.E_ionizations
    end

    AURORA.load_or_compute_cascading!(AURORA.CascadingCache(), energy_grid; policy=load_policy)
    @test compatible_cache_count(n2_dir) >= 1

    AURORA.load_or_compute_cascading!(AURORA.CascadingCache(), energy_grid; policy=skip_save_policy)
    skip_n2_dir = joinpath(cache_root, "skip_save", "e_cascading", "N2")
    @test isempty(cache_files(skip_n2_dir))

    AURORA.clear_cascading_cache!(cache_root=cache_root)
    @test isempty(cache_files(n2_dir))
    @test isempty(cache_files(o2_dir))
    @test isempty(cache_files(o_dir))
end
