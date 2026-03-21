using Dates: Dates, now
using MAT: matopen
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

function ScatteringData(θ_lims; n_direction=720, verbose=true)
    validate_θ_lims(θ_lims)
    Ω_beam = beam_weight(θ_lims)
    P_scatter, Ω_subbeam_relative, θ1 = find_scattering_matrices(θ_lims, n_direction;
                                                                 verbose)
    FT = eltype(P_scatter)
    return ScatteringData{FT, typeof(P_scatter), typeof(Ω_subbeam_relative), typeof(Ω_beam)}(
        P_scatter, Ω_subbeam_relative, Ω_beam, vec(θ1)
    )
end

"""
    load_scattering_matrices(θ_lims)

Load the scattering matrices for the given pitch-angle limits.

# Calling
`μ_lims, μ_center, scattering = load_scattering_matrices(θ_lims)`

# Inputs
- `θ_lims`: pitch angle limits of the e- beams (deg). Vector [n_beam + 1]

# Outputs
- `μ_lims`: cosine of the pitch angle limits of the e- beams. Vector [n_beam + 1]
- `μ_center`: cosine of the pitch angle of the middle of the e- beams. Vector [n_beam]
- `scattering`: Tuple with several of the scattering informations, namely scattering
    = `(P_scatter, Ω_subbeam_relative, Ω_beam)`
    + `P_scatter`: probabilities for scattering in 3D from beam to beam. Matrix [n`_`direction x n`_`direction]
    + `Ω_subbeam_relative`: relative contribution from within each beam. Matrix [n`_`beam x n`_`direction]
    + `Ω_beam`: solid angle for each stream (ster). Vector [n_beam]
    + `θ_scatter`: scattering angles used in the calculations. Vector [n_direction]
"""
function load_scattering_matrices(θ_lims)
    validate_θ_lims(θ_lims)
    μ_lims = cosd.(θ_lims)
    μ_center = mu_avg(θ_lims)
    Ω_beam = beam_weight(θ_lims) # this beam weight is calculated in a continuous way
    P_scatter, Ω_subbeam_relative, θ₁ = find_scattering_matrices(θ_lims, 720)
    scattering = (P_scatter = P_scatter, Ω_subbeam_relative = Ω_subbeam_relative,
                     Ω_beam = Ω_beam, θ_scatter = θ₁)

    return μ_lims, μ_center, scattering
end

"""
    validate_θ_lims(θ_lims)

Validate that the pitch-angle limits `θ_lims` are correctly specified.
Throws an `ArgumentError` if:
- `θ_lims` does not include 180° (field-aligned downward)
- `θ_lims` does not include 0° (field-aligned upward)
- `θ_lims` is not in descending order

# Inputs
- `θ_lims`: pitch-angle limits of the electron beams (e.g. 180:-10:0)
"""
function validate_θ_lims(θ_lims)
    if maximum(θ_lims) != 180
        throw(ArgumentError("θ_lims must include 180° (field-aligned downward). " *
              "Got maximum of $(maximum(θ_lims))°. " *
              "Example of valid input: 180:-10:0."))
    end
    if minimum(θ_lims) != 0
        throw(ArgumentError("θ_lims must include 0° (field-aligned upward). " *
              "Got minimum of $(minimum(θ_lims))°. " *
              "Example of valid input: 180:-10:0."))
    end
    if !issorted(θ_lims, rev=true)
        throw(ArgumentError("θ_lims must be in descending order (e.g., 180:-10:0). " *
              "Got: $θ_lims"))
    end
    return nothing
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
    find_scattering_matrices(θ_lims, n_direction=720)

Look for scattering matrices that match the pitch-angle limits `θ_lims` and the number
of direction/sub-beams `n_direction`. If a file is found, the scattering matrices are
directly loaded. Otherwise, they are calculated and saved to a file.

# Calling
`P_scatter, Ω_subbeam_relative, θ₁ = find_scattering_matrices(θ_lims, n_direction)`

# Inputs
- `θ_lims`: pitch-angle limits of the electron beams (e.g. 180:-10:0), where 180°
    corresponds to field aligned down, and 0° field aligned up.
- `n_direction`: number of directions or sub-beams to use for the discretized calculations
    of the scattering matrices. Defaults to 720 when left empty.

# Outputs
- `P_scatter`: probabilities for scattering in 3D from beam to beam. Matrix [n`_`direction x
    n`_`direction]
- `Ω_subbeam_relative`: relative weight of each sub-beam within each beam, normalized so that
    summing along the sub-beams gives 1 for each beam. Matrix [n`_`beam x n`_`direction]
- `θ₁`: scattering angles used in the calculations. Vector [n_direction]
"""
function find_scattering_matrices(θ_lims, n_direction=720; verbose = true)
    scattering_files = readdir(pkgdir(AURORA, "internal_data", "e_scattering"))
    found_them = 0
    for i1 in eachindex(scattering_files)
        if !isdir(scattering_files[i1])
            try
                filename = pkgdir(AURORA, "internal_data", "e_scattering", scattering_files[i1])
                file = matopen(filename)
                θ_lims_file = read(file, "theta_lims")
                n_direction_file = read(file, "n_direction")
                close(file)
                if θ_lims == θ_lims_file && n_direction == n_direction_file
                    file = matopen(filename)
                    P_scatter = read(file, "P_scatter")
                    Ω_subbeam_relative = read(file, "subbeamweight_relative")
                    θ₁ = read(file, "theta_scatter")
                    close(file)
                    found_them = 1
                    verbose && println("Loading scattering-matrices from file: ", scattering_files[i1], " ✅")
                    break
                end
            catch
            end
        end
    end
    if found_them == 0
        verbose && println("Could not find file with matching pitch-angle grid.")
        verbose && println("Starting to calculate the requested scattering-matrices.")

        P_scatter, Ω_subbeam_relative, θ₁ = calculate_scattering_matrices(θ_lims, n_direction)

        # Save the results for future use using the current internal naming scheme.
        filename = pkgdir(AURORA, "internal_data", "e_scattering",
                            string(length(θ_lims) - 1, "_streams_",
                            Dates.format(now(), "yyyymmdd-HHMMSS"),
                            ".mat"))
        file = matopen(filename, "w")
        write(file, "P_scatter", P_scatter)
        write(file, "subbeamweight_relative", Ω_subbeam_relative)
        write(file, "theta_scatter", θ₁)
        write(file, "theta_lims", Vector(θ_lims))
        write(file, "n_direction", n_direction)
        close(file)
    end

    return P_scatter, Ω_subbeam_relative, θ₁
end






# Another 100x factor (0.2s vs 20s for the legacy version, still with θ_lims = 180:-10:0 and
# n_direction = 720
"""
    calculate_scattering_matrices(θ_lims, n_direction = 720)

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
function calculate_scattering_matrices(θ_lims, n_direction = 720)
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
    p = Progress(length(θ₀); desc=string("Calculating scattering matrices "), color=:blue);
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
        next!(p)
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
    calculate_scattering_matrices_legacy(θ_lims, n_direction=720)

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
function calculate_scattering_matrices_legacy(θ_lims, n_direction=720)
    μ_lims = cosd.(θ_lims)

    θ₀ = Vector(0:(180/n_direction):180); θ₀ = deg2rad.((θ₀[1:end - 1] .+ θ₀[2:end]) / 2)
    θ₁ = deg2rad.(0:(180/n_direction):360); θ₁ = θ₁[1:end-1]
    ϕ = deg2rad.(0:(180/n_direction):360); ϕ = ϕ[1:end-1]

    μ = Vector{Float64}(undef, length(ϕ))
    B = Array{Float64}(undef, length(θ₀), length(θ₁), length(θ_lims) - 1)

    st = sin.(ϕ)
    ct = cos.(ϕ)
    cache = Vector{Float64}(undef, length(ϕ))

    p = Progress(length(θ₀); desc=string("Calculating scattering matrices "), color=:blue);

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
        next!(p)
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
