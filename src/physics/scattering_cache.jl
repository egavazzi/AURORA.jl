using Dates: Dates, now
using JLD2: jldopen


"""
    load_or_compute_scattering_cache(θ_lims, n_direction=720;
                                     verbose=true, policy=CachePolicy())

Look for scattering matrices that match the pitch-angle limits `θ_lims` and the number
of direction/sub-beams `n_direction`. If a cache file is found, the scattering matrices are
directly loaded. Otherwise, they are calculated and optionally saved to a file.

# Calling
`P_scatter, Ω_subbeam_relative, θ₁ = load_or_compute_scattering_cache(θ_lims, n_direction)`

# Inputs
- `θ_lims`: pitch-angle limits of the electron beams (e.g. 180:-10:0), where 180°
    corresponds to field aligned down, and 0° field aligned up.
- `n_direction`: number of directions or sub-beams to use for the discretized calculations
    of the scattering matrices. Defaults to 720 when left empty.
- `policy.cache_root`: parent directory that contains the `e_scattering/` cache subtree

# Outputs
- `P_scatter`: probabilities for scattering in 3D from beam to beam. Matrix [n`_`direction x
    n`_`direction]
- `Ω_subbeam_relative`: relative weight of each sub-beam within each beam, normalized so that
    summing along the sub-beams gives 1 for each beam. Matrix [n`_`beam x n`_`direction]
- `θ₁`: scattering angles used in the calculations. Vector [n_direction]
"""
function load_or_compute_scattering_cache(θ_lims, n_direction=720;
                                          verbose = true,
                                          policy::CachePolicy = CachePolicy())
    file_found, filepath = policy.force_recompute ? (false, "") :
                           find_scattering_cache_file(θ_lims, n_direction; verbose, policy)

    if file_found
        try
            return load_scattering_cache(filepath; verbose)
        catch err
            @warn "Failed to load scattering cache $(basename(filepath)). Recomputing." exception = err
        end
    end

    if !file_found && !policy.force_recompute
        verbose && println("No compatible scattering cache found. Computing...")
    end

    P_scatter, Ω_subbeam_relative, θ₁ = calculate_scattering_matrices(θ_lims, n_direction; verbose)

    if policy.save_cache
        save_scattering_cache(P_scatter, Ω_subbeam_relative, θ₁, θ_lims, n_direction; verbose, policy)
    else
        verbose && println("Scattering cache not saved (save_cache=false).")
    end

    return P_scatter, Ω_subbeam_relative, θ₁
end

function find_scattering_cache_file(θ_lims, n_direction;
                                    verbose::Bool = true,
                                    policy::CachePolicy = CachePolicy())
    cache_dir = scattering_cache_dir(policy)
    isdir(cache_dir) || return (false, "")
    for filename in readdir(cache_dir)
        endswith(filename, ".jld2") || continue

        filepath = joinpath(cache_dir, filename)
        isdir(filepath) && continue

        try
            result = jldopen(filepath, "r") do file
                required_keys = ["version_AURORA", "P_scatter", "Ω_subbeam_relative",
                                 "theta_scatter", "theta_lims", "n_direction"]
                all(haskey(file, key) for key in required_keys) || return nothing

                version_saved = string(file["version_AURORA"])
                if version_saved != cache_version_string()
                    verbose && println("Skipping $(filename): built with AURORA $version_saved.")
                    return nothing
                end

                θ_lims_file = file["theta_lims"]
                n_direction_file = file["n_direction"]
                if θ_lims == θ_lims_file && n_direction == n_direction_file
                    return (true, filepath)
                end
                return nothing
            end
            isnothing(result) || return result
        catch
            continue
        end
    end

    return (false, "")
end

function load_scattering_cache(filepath; verbose = true)
    verbose && println("Loading scattering cache from file: $(basename(filepath))")

    P_scatter = nothing
    Ω_subbeam_relative = nothing
    θ_scatter = nothing
    jldopen(filepath, "r") do file
        version_saved = string(file["version_AURORA"])
        if version_saved != cache_version_string()
            error("Found incompatible scattering cache file $(basename(filepath)). Recomputing.")
        end
        P_scatter = file["P_scatter"]
        Ω_subbeam_relative = file["Ω_subbeam_relative"]
        θ_scatter = file["theta_scatter"]
    end

    return P_scatter, Ω_subbeam_relative, θ_scatter
end

function save_scattering_cache(P_scatter, Ω_subbeam_relative, θ_scatter, θ_lims, n_direction;
                               verbose::Bool = true,
                               policy::CachePolicy = CachePolicy())
    cache_dir = scattering_cache_dir(policy)
    mkpath(cache_dir)
    filename = joinpath(cache_dir,
                        string("scattering_", length(θ_lims) - 1, "_streams_",
                               Dates.format(now(), "yyyymmdd-HHMMSS"), ".jld2"))
    jldopen(filename, "w") do file
        file["version_AURORA"] = cache_version_string()
        file["P_scatter"] = P_scatter
        file["Ω_subbeam_relative"] = Ω_subbeam_relative
        file["theta_scatter"] = θ_scatter
        file["theta_lims"] = Vector(θ_lims)
        file["n_direction"] = n_direction
    end
    verbose && println("Saved scattering cache to $(basename(filename)).")
    return filename
end

function clear_scattering_cache!(; cache_root::String = default_cache_root())
    cache_dir = joinpath(cache_root, "e_scattering")
    isdir(cache_dir) || return nothing

    for filename in readdir(cache_dir)
        endswith(filename, ".jld2") || continue
        rm(joinpath(cache_dir, filename); force = true)
    end

    return nothing
end

function scattering_cache_dir(policy::CachePolicy = CachePolicy())
    return joinpath(policy.cache_root, "e_scattering")
end
