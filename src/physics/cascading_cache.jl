using Dates: Dates, now
using JLD2: jldopen

# Keys a cascading cache file must carry. Besides the matrices themselves, the file records
# every ingredient of the physics it was built from, so that a file built for a different
# spec (other thresholds, other secondary counts, other secondary law) is rejected instead
# of being reused under the same species name.
const CASCADING_CACHE_KEYS = ("version_AURORA", "Q_primary", "Q_secondary", "E_edges",
                              "E_ionizations", "n_secondaries", "law_fingerprint")

"""
    load_or_compute_cascading!(cache::SpeciesCascadingCache, energy_grid; verbose=true, policy=CachePolicy())

Populate the cascading transfer matrices for one neutral species.

Matrices are loaded from a JLD2 cache file when one on disk was built with the same
`version_AURORA`, a grid whose edges start with the requested ones, and the same
ionization thresholds, secondary counts and secondary law as `cache.spec`. Otherwise they
are computed from scratch (and optionally saved to disk depending on `CachePolicy` options).

A spec whose secondary law is not [`is_fingerprintable`](@ref) is never cached: there is
nothing to tell one such law from another, so its matrices are recomputed every time.
"""
function load_or_compute_cascading!(cache::SpeciesCascadingCache, energy_grid::EnergyGrid;
                                    verbose::Bool = true,
                                    policy::CachePolicy = CachePolicy())
    E_edges = energy_grid.E_edges
    n_E = length(E_edges)
    cacheable = is_fingerprintable(cache.spec.secondary_law)

    file_found, filepath = (policy.force_recompute || !cacheable) ? (false, "") :
                           find_cascading_cache(cache.spec, E_edges; verbose, policy)

    if file_found
        try
            cascading_data = load_cascading_cache(filepath, n_E, cache.spec; verbose)
            cache.primary_transfer_matrix = cascading_data[1]
            cache.secondary_transfer_matrix = cascading_data[2]
            cache.E_edges = cascading_data[3]
            # The file's thresholds were checked for equality with the spec's before it was
            # accepted, so the spec is the single source of truth for them.
            cache.ionization_thresholds = copy(cache.spec.ionization_thresholds)
            return nothing
        catch err
            @warn "Could not load cascading cache $(basename(filepath)). Recomputing." exception = err
        end
    end

    if !cacheable
        verbose && println("Cascading matrices for $(cache.spec.name) are not cached: its \
                            secondary law is a $(typeof(cache.spec.secondary_law)), which \
                            carries nothing to identify it by. Wrap the law with @law, or \
                            hold its parameters in a functor struct, to cache them. \
                            Computing...")
    elseif !file_found && !policy.force_recompute
        verbose && println("No compatible cascading cache for $(cache.spec.name). Computing...")
    end

    cascading_data = calculate_cascading_matrices(cache.spec, E_edges; verbose)
    cache.primary_transfer_matrix = cascading_data[1]
    cache.secondary_transfer_matrix = cascading_data[2]
    cache.E_edges = cascading_data[3]
    cache.ionization_thresholds = cascading_data[4]

    if cacheable
        if policy.save_cache
            save_cascading_cache(cache; verbose, policy)
        else
            verbose && println("Cascading cache for $(cache.spec.name) not saved (save_cache=false).")
        end
    end

    return nothing
end

function find_cascading_cache(spec::CascadingSpec, E_edges;
                              verbose::Bool = true,
                              policy::CachePolicy = CachePolicy())
    species_dir = cascading_cache_dir(spec, policy)
    isdir(species_dir) || return (false, "")
    cascading_files = readdir(species_dir)
    fingerprint = law_fingerprint(spec.secondary_law)

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
                if !all(haskey(file, key) for key in CASCADING_CACHE_KEYS)
                    verbose && println("Skipping $(filename): missing cache entries.")
                    return nothing
                end

                version_saved = file["version_AURORA"]
                if string(version_saved) != cache_version_string()
                    verbose && println("Skipping $(filename): built with AURORA $version_saved.")
                    return nothing
                end

                if file["E_ionizations"] != spec.ionization_thresholds
                    verbose && println("Skipping $(filename): built for ionization thresholds \
                                        $(file["E_ionizations"]).")
                    return nothing
                end

                if file["n_secondaries"] != spec.n_secondaries
                    verbose && println("Skipping $(filename): built for secondary counts \
                                        $(file["n_secondaries"]).")
                    return nothing
                end

                if file["law_fingerprint"] != fingerprint
                    verbose && println("Skipping $(filename): built with a different \
                                        secondary law.")
                    return nothing
                end

                E_edges_saved = file["E_edges"]
                if length(E_edges) <= length(E_edges_saved) && E_edges_saved[1:length(E_edges)] == E_edges
                    return (true, filepath)
                end
            end
            isnothing(result) || return result
        catch err
            verbose && println("Skipping $(filename): $(err)")
            continue
        end
    end

    return (false, "")
end

function load_cascading_cache(filepath, n_E::Int, spec::CascadingSpec; verbose::Bool = true)
    verbose && println("Loading cascading matrices from file: $(basename(filepath))")
    n_bins = n_E - 1
    n_thresholds = length(spec.ionization_thresholds)
    return jldopen(filepath, "r") do file
        Q_primary     = file["Q_primary"][1:n_bins, 1:n_bins, :]
        Q_secondary   = file["Q_secondary"][1:n_bins, 1:n_bins, :]
        E_edges       = file["E_edges"][1:n_E]
        E_ionizations = file["E_ionizations"]
        if size(Q_primary, 3) != n_thresholds || size(Q_secondary, 3) != n_thresholds
            throw(ArgumentError(
                "cascading cache $(basename(filepath)) holds $(size(Q_primary, 3)) \
                 degraded-primary and $(size(Q_secondary, 3)) secondary matrices, but \
                 $(spec.name) has $(n_thresholds) ionization thresholds"))
        end
        (Q_primary, Q_secondary, E_edges, E_ionizations)
    end
end

function save_cascading_cache(cache::SpeciesCascadingCache;
                              verbose::Bool = true,
                              policy::CachePolicy = CachePolicy())
    species_dir = cascading_cache_dir(cache.spec, policy)
    mkpath(species_dir)
    filename = joinpath(species_dir,
                       string("cascading_", cache.spec.name, "_",
                             Dates.format(now(), "yyyymmdd-HHMMSS"),
                             ".jld2"))
    jldopen(filename, "w") do file
        file["version_AURORA"]  = cache_version_string()
        file["Q_primary"]       = cache.primary_transfer_matrix
        file["Q_secondary"]     = cache.secondary_transfer_matrix
        file["E_edges"]         = cache.E_edges
        file["E_ionizations"]   = cache.ionization_thresholds
        file["n_secondaries"]   = cache.spec.n_secondaries
        file["law_fingerprint"] = law_fingerprint(cache.spec.secondary_law)
    end
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
