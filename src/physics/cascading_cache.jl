using Dates: Dates, now
using JLD2: jldopen, @load, @save


"""
    load_or_compute_cascading_cache!(cache, energy_grid; verbose=true, policy=CachePolicy())

Populate the cascading cache for all neutral species.

If compatible (i.e. same grid, same `version_AURORA`, etc) JLD2 cache files are found on
disk, the matrices are loaded from there. If not, they are computed from scratch (and
possibly saved to disk depending on `CachePolicy` options).
"""
function load_or_compute_cascading_cache!(cache::CascadingCache, energy_grid::EnergyGrid;
                                          verbose::Bool = true,
                                          policy::CachePolicy = CachePolicy())
    for species_cache in cache
        load_or_compute_cascading_cache!(species_cache, energy_grid; verbose, policy)
    end
    return nothing
end

function load_or_compute_cascading_cache!(cache::SpeciesCascadingCache, energy_grid::EnergyGrid;
                                          verbose::Bool = true,
                                          policy::CachePolicy = CachePolicy())
    E_edges = energy_grid.E_edges

    file_found, filepath = policy.force_recompute ? (false, "") :
                           find_cascading_cache_file(cache.spec, E_edges; verbose, policy)

    if file_found
        try
            cascading_data = load_cascading_cache(filepath; verbose)
            cache.primary_transfer_matrix = cascading_data[1][1:length(E_edges)-1, 1:length(E_edges)-1, :]
            cache.secondary_transfer_matrix = cascading_data[2][1:length(E_edges)-1, 1:length(E_edges)-1, :]
            cache.E_edges = cascading_data[3][1:length(E_edges)]
            cache.ionization_thresholds = cascading_data[4]
            return nothing
        catch err
            @warn "Failed to load cascading cache $(basename(filepath)). Recomputing." exception = err
        end
    end

    if !file_found && !policy.force_recompute
        verbose && println("No compatible cascading cache for $(cache.spec.name). Computing...")
    end

    cascading_data = calculate_cascading_matrices(cache.spec, E_edges; verbose)
    cache.primary_transfer_matrix = cascading_data[1]
    cache.secondary_transfer_matrix = cascading_data[2]
    cache.E_edges = cascading_data[3]
    cache.ionization_thresholds = cascading_data[4]

    if policy.save_cache
        save_cascading_cache(cache.primary_transfer_matrix, cache.secondary_transfer_matrix,
                             cache.E_edges, cache.ionization_thresholds,
                             cache.spec.name; verbose, policy)
    else
        verbose && println("Cascading cache for $(cache.spec.name) not saved (save_cache=false).")
    end

    return nothing
end

function find_cascading_cache_file(spec::CascadingSpec, E_edges;
                                   verbose::Bool = true,
                                   policy::CachePolicy = CachePolicy())
    species_dir = cascading_cache_dir(spec, policy)
    isdir(species_dir) || return (false, "")
    cascading_files = readdir(species_dir)

    for filename in cascading_files
        if !endswith(filename, ".jld2")
            continue
        end

        if isdir(joinpath(species_dir, filename))
            continue
        end

        filepath = joinpath(species_dir, filename)
        try
            result = jldopen(filepath, "r") do file
                required_keys = ["version_AURORA", "Q_primary", "Q_secondary", "E_edges", "E_ionizations"]
                if !all(haskey(file, key) for key in required_keys)
                    return nothing
                end

                version_saved = file["version_AURORA"]
                if string(version_saved) != cache_version_string()
                    verbose && println("Skipping $(filename): built with AURORA $version_saved.")
                    return nothing
                end

                E_edges_saved = file["E_edges"]
                if length(E_edges) <= length(E_edges_saved) && E_edges_saved[1:length(E_edges)] == E_edges
                    return (true, filepath)
                end
            end
            isnothing(result) || return result
        catch
            continue
        end
    end

    return (false, "")
end

function load_cascading_cache(filepath; verbose::Bool = true)
    verbose && println("Loading cascading matrices from file: $(basename(filepath))")
    @load filepath Q_primary Q_secondary E_ionizations E_edges
    return (Q_primary, Q_secondary, E_edges, E_ionizations)
end

function save_cascading_cache(primary_transfer_matrix, secondary_transfer_matrix,
                              E_edges, ionization_thresholds, species_name;
                              verbose::Bool = true,
                              policy::CachePolicy = CachePolicy())
    species_dir = cascading_cache_dir(CascadingSpec(species_name, Float64[], (_, _) -> 0.0),
                                      policy)
    mkpath(species_dir)
    filename = joinpath(species_dir,
                       string("cascading_", species_name, "_",
                             Dates.format(now(), "yyyymmdd-HHMMSS"),
                             ".jld2"))
    version_AURORA = cache_version_string()
    Q_primary = primary_transfer_matrix
    Q_secondary = secondary_transfer_matrix
    E_ionizations = ionization_thresholds
    @save filename version_AURORA Q_primary Q_secondary E_edges E_ionizations
    verbose && println("Saved cascading matrices to $(basename(filename)).")
    return filename
end

function clear_cascading_cache!(; cache_root::String = default_cache_root())
    base_dir = joinpath(cache_root, "e_cascading")
    isdir(base_dir) || return nothing
    for species_name in readdir(base_dir)
        species_dir = joinpath(base_dir, species_name)
        isdir(species_dir) || continue
        for filename in readdir(species_dir)
            endswith(filename, ".jld2") || continue
            rm(joinpath(species_dir, filename); force = true)
        end
    end
    return nothing
end

function cascading_cache_dir(spec::CascadingSpec, policy::CachePolicy = CachePolicy())
    return joinpath(policy.cache_root, "e_cascading", spec.name)
end
