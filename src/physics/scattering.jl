using Dates: Dates, now
using JLD2: jldopen
using ProgressMeter: Progress, next!
using Rotations: AngleAxis


"""
    ScatteringData{FT, A<:AbstractArray{FT}, M<:AbstractMatrix{FT}, V<:AbstractVector{FT}}

Pre-computed scattering probability matrices for beam-to-beam scattering.
"""
struct ScatteringData{FT, A<:AbstractArray{FT}, M<:AbstractMatrix{FT}, V<:AbstractVector{FT}}
    P_scatter::A
    Ω_subbeam_relative::M
    Ω_beam::V
    θ_scatter::V
end

function ScatteringData(grid::PitchAngleGrid;
                        n_direction = 720, verbose = true,
                        policy::CachePolicy = CachePolicy())
    θ_lims = grid.θ_lims
    Ω_beam = beam_weight(θ_lims)
    P_scatter, Ω_subbeam_relative, θ1 = load_or_compute_scattering_cache(θ_lims, n_direction;
                                                                         verbose, policy)
    FT = eltype(P_scatter)
    return ScatteringData{FT, typeof(P_scatter), typeof(Ω_subbeam_relative), typeof(Ω_beam)}(
        P_scatter, Ω_subbeam_relative, Ω_beam, vec(θ1)
    )
end

function Base.show(io::IO, sd::ScatteringData)
    n_beams = length(sd.Ω_beam)
    n_dir = length(sd.θ_scatter)
    print(io, "ScatteringData($(n_beams) beams, $(n_dir) directions)")
end

function Base.show(io::IO, ::MIME"text/plain", sd::ScatteringData)
    n_beams = length(sd.Ω_beam)
    n_dir = length(sd.θ_scatter)
    println(io, "ScatteringData:")
    println(io, "├── Beams: $(n_beams)")
    println(io, "├── Directions: $(n_dir)")
    print(io, "└── P_scatter: $(join(size(sd.P_scatter), 'x'))")
end


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
function scattering_cache_dir(policy::CachePolicy = CachePolicy())
    return joinpath(policy.cache_root, "e_scattering")
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
            error("Found incompatible scattering cache file $(basename(filepath)); recomputing.")
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

"""
    clear_scattering_cache!(; cache_root=default_cache_root())

Remove generated JLD2 scattering cache files from the scattering cache directory under
`cache_root`.
"""
function clear_scattering_cache!(; cache_root::String = default_cache_root())
    cache_dir = joinpath(cache_root, "e_scattering")
    isdir(cache_dir) || return nothing

    for filename in readdir(cache_dir)
        endswith(filename, ".jld2") || continue
        rm(joinpath(cache_dir, filename); force = true)
    end

    return nothing
end

function load_or_compute_scattering_cache(θ_lims, n_direction=720;
                                          verbose = true,
                                          policy::CachePolicy = CachePolicy())
    file_found, filepath = policy.force_recompute ? (false, "") :
                           find_scattering_cache_file(θ_lims, n_direction; verbose, policy)

    if file_found
        try
            return load_scattering_cache(filepath; verbose)
        catch err
            @warn "Failed to load scattering cache $(basename(filepath)); recomputing."
        end
    end

    if !file_found && !policy.force_recompute
        verbose && println("No compatible scattering cache found; computing...")
    end

    P_scatter, Ω_subbeam_relative, θ₁ = calculate_scattering_matrices(θ_lims, n_direction; verbose)

    if policy.save_cache
        save_scattering_cache(P_scatter, Ω_subbeam_relative, θ₁, θ_lims, n_direction; verbose, policy)
    else
        verbose && println("Scattering cache not saved (save_cache=false).")
    end

    return P_scatter, Ω_subbeam_relative, θ₁
end






# Another 100x factor (0.2s vs 20s for the legacy version, still with θ_lims = 180:-10:0 and
# n_direction = 720
"""
    calculate_scattering_matrices(θ_lims, n_direction = 720; verbose = true)

Calculate the scattering matrices for given pitch-angle limits `θ_lims` of the electron
beams. Uses 720 directions by default for the start angle and the scattering angles,
equivalent to (1/4)° steps.

# Calling
`P_scatter, Ω_subbeam_relative, θ₁ = calculate_scattering_matrices(θ_lims, n_direction)`

# Inputs
- `θ_lims`: pitch-angle limits of the electron beams (e.g. 180:-10:0), where 180°
    corresponds to field aligned down, and 0° field aligned up.
- `n_direction`: number of directions or sub-beams to use for the discretized calculations
    of the scattering matrices. Defaults to 720 when left empty.

# Outputs
- `P_scatter`: probabilities for scattering in 3D from beam to beam. Matrix [n\\_direction x n\\_direction x n\\_beam]
- `Ω_subbeam_relative`: relative weight of each sub-beam within each beam, normalized so that
    summing along the sub-beams gives 1 for each beam. Matrix [n\\_beam x n\\_direction]
- `θ₁`: scattering angles used in the calculations. Vector [n\\_direction]
"""
function calculate_scattering_matrices(θ_lims, n_direction = 720; verbose = true)
    μ_lims = cosd.(θ_lims)

    θ₀ = Vector(0:(180/n_direction):180); θ₀ = deg2rad.((θ₀[1:end - 1] .+ θ₀[2:end]) / 2)
    θ₁ = deg2rad.(0:(180/n_direction):180); θ₁ = θ₁[1:end-1]
    ct₀ = cos.(θ₀)
    st₀ = sin.(θ₀)
    ct₁ = cos.(θ₁)
    st₁ = sin.(θ₁)

    P_scatter = Array{Float64}(undef, length(θ₀), length(θ₁), length(θ_lims) - 1)
    Ω_subbeam = Array{Float64}(undef, length(θ_lims) - 1, length(θ₀))

    #=
    What this does:
    - Take an initial angle θ₀
    - Take a scattering angle θ₁
    - Instead of manually doing a rotation around the vector defined by θ₀ and then taking
        the z-component of each point, calculate analytically the distribution of the
        projection on the z-axis of the scattering "circle". This is much faster (and more
        accurate?)
    =#
    p = verbose ? Progress(length(θ₀); desc="Calculating scattering matrices ", color=:blue) : nothing
    for i0 in length(θ₀):-1:1
        for i1 in length(θ₁):-1:1
            for iμ in (length(μ_lims) - 1):-1:1
                # This is an analytical calculation of the distribution of the projection
                # of the scattered vectors on the z-axis.
                z₁ = μ_lims[iμ]
                z₂ = μ_lims[iμ + 1]
                P_scatter[i0, i1, iμ] = 1 / π *
                                  (asin(clamp((z₂ - ct₀[i0]ct₁[i1]) / (st₀[i0]st₁[i1]), -1, 1)) -
                                   asin(clamp((z₁ - ct₀[i0]ct₁[i1]) / (st₀[i0]st₁[i1]), -1, 1)))

            end
        end
        verbose && next!(p)
    end

    # Normalize so that all sum(P_scatter[i, j, :], dims=3) = 1
    P_scatter = P_scatter ./ repeat(sum(P_scatter, dims = 3), outer = (1, 1, size(P_scatter, 3)))
    # Calculate the beam weight of each subdivision
    for iμ in (length(μ_lims) - 1):-1:1
        Ω_subbeam[iμ, :] .= abs.(sin.(θ₀)) .* (μ_lims[iμ] .< cos.(θ₀) .< μ_lims[iμ + 1])
    end
    # Normalize so that all sum(Ω_subbeam_relative, dims=2) = 1
    beam_sums = sum(Ω_subbeam, dims=2)
    if any(iszero, beam_sums)
        error(
            "Some pitch-angle beams have no sub-beam coverage. " *
            "This typically happens when θ_lims doesn't span the full range from 180° to 0°. " *
            "Make sure θ_lims includes both 180° and 0°."
        )
    end
    Ω_subbeam_relative = Ω_subbeam ./ repeat(beam_sums, 1, size(Ω_subbeam, 2))

    return P_scatter, Ω_subbeam_relative, θ₁
end












############################################################################################
##                                       Legacy                                           ##
############################################################################################

# 144x faster than the Matlab code (24s instead of 1h for θ_lims = 180:-10:0 and
# n_direction = 720)
"""
    calculate_scattering_matrices_legacy(θ_lims, n_direction=720; verbose = true)

Calculate the scattering matrices for given pitch-angle limits `θ_lims` of the electron
beams.

# Calling
`Pmu2mup, theta2beamW, BeamWeight_relative, θ₁ =  calculate_scattering_matrices_legacy(θ_lims, n_direction)`

# Inputs
- `θ_lims`: pitch-angle limits of the electron beams (e.g. 180:-10:0), where 180°
    corresponds to field aligned down, and 0° field aligned up.
- `n_direction`: number of directions or sub-beams to use for the discretized calculations
    of the scattering matrices. Defaults to 720 when left empty.

# Outputs
- `Pmu2mup`: probabilities for scattering in 3D from beam to beam. Matrix [n\\_direction x
    2 * n\\_direction x n\\_beam]
- `theta2beamW`: weight of each sub-beam within each beam. Matrix [n\\_beam x
    n\\_direction]
- `BeamWeight_relative`: relative weight of each sub-beam within each beam. It is the same
    as theta2beamW but normalized so that summing along the sub-beams gives 1 for each beam.
    Matrix [n\\_beam x n\\_direction]
- `θ₁`: scattering angles used in the calculations. Vector [n_direction]
"""
function calculate_scattering_matrices_legacy(θ_lims, n_direction=720; verbose = true)
    μ_lims = cosd.(θ_lims)

    θ₀ = Vector(0:(180/n_direction):180); θ₀ = deg2rad.((θ₀[1:end - 1] .+ θ₀[2:end]) / 2)
    θ₁ = deg2rad.(0:(180/n_direction):360); θ₁ = θ₁[1:end-1]
    ϕ = deg2rad.(0:(180/n_direction):360); ϕ = ϕ[1:end-1]

    μ = Vector{Float64}(undef, length(ϕ))
    B = Array{Float64}(undef, length(θ₀), length(θ₁), length(θ_lims) - 1)

    st = sin.(ϕ)
    ct = cos.(ϕ)
    cache = Vector{Float64}(undef, length(ϕ))

    p = verbose ? Progress(length(θ₀); desc="Calculating scattering matrices ", color=:blue) : nothing

    rotation_axis1 = [0, 1, 0]
    #=
    What this basically do is:
        - start with an initial angle θ₀ (from 0 to 180°), give us a vector e0
        - rotate one time by a scattering angle θ₁ (from 0 to 360°), give us a vector e1
        - rotate the vector e1 around e0 by an angle ϕ from 0 to 360°, each time giving us a vector es
        - save the "z"-position of the es vector, which gives us the final "beam"-position
    =#
    for i0 in length(θ₀):-1:1
        e0 = [sin(θ₀[i0]), 0, cos(θ₀[i0])]

        for i1 in length(θ₁):-1:1
            rotation_angle = θ₁[i1]
            rotation_matrix = AngleAxis(rotation_angle, rotation_axis1[1], rotation_axis1[2], rotation_axis1[3])
            e1 = rotation_matrix * e0

            for iϕ in length(ϕ):-1:1
                rotation_axis = e0
                μ[iϕ] = my_very_special_rotation(e1, ct[iϕ], st[iϕ], rotation_axis)
                # rotation_angle = ϕ[iϕ]
                # rotation_matrix = AngleAxis(rotation_angle, rotation_axis[1], rotation_axis[2], rotation_axis[3])
                # es = rotation_matrix * e1
                # μ[iϕ] = es[3]
            end

            for iμ in (length(μ_lims) - 1):-1:1
                cache .= μ_lims[iμ] .< μ .< μ_lims[iμ + 1]
                B[i0, i1, iμ] = sum(cache)
            end
        end
        verbose && next!(p)
    end

    Pmu2mup = B ./ repeat(sum(B, dims=3), outer = (1, 1, size(B, 3)))

    theta2beamW = Array{Float64}(undef, length(θ_lims) - 1, length(θ₀))
    for iμ in (length(μ_lims) - 1):-1:1
        theta2beamW[iμ, :] .= abs.(sin.(θ₀)) .* (μ_lims[iμ] .< cos.(θ₀) .< μ_lims[iμ + 1])
    end
    BeamWeight_relative = theta2beamW ./ repeat(sum(theta2beamW, dims=2), 1, size(theta2beamW, 2))

    return Pmu2mup, theta2beamW, BeamWeight_relative, θ₁
end

function my_very_special_rotation(v, st, ct, k)
    # Rotate the vector`v` around the rotation_axis `k` by a rotation_angle ϕ defined by
    # cos(ϕ) = `ct` and sin(ϕ) = `st`. It is special in the way that it only computes and
    # returns the z-component of the final vector.
    # See Rodrigues' rotation formula at
    # https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula#Matrix_notation
    return v[3] * ct + (k[1]v[2] - k[2]v[1]) * st + k[3] * (k[1]v[1] + k[2]v[2] + k[3]v[3]) * (1 - ct)
end
