using Dates: Dates, now
using MAT: matopen
using ProgressMeter: Progress, next!
using Rotations: AngleAxis

"""
    find_scattering_matrices(θ_lims, n_direction=720)

Look for scattering matrices that match the pitch-angle limits `θ_lims` and the number
of direction/sub-beams `n_direction`. If a file is found, the scattering matrices are
directly loaded. Otherwise, they are calculated and saved to a file.

# Calling
`Pmu2mup, theta2beamW, BeamWeight_relative, θ₁ =  find_scattering_matrices(θ_lims, n_direction)`

# Inputs
- `θ_lims`: pitch-angle limits of the electron beams (e.g. 180:-10:0), where 180°
    corresponds to field aligned down, and 0° field aligned up.
- `n_direction`: number of directions or sub-beams to use for the discretized calculations
    of the scattering matrices. Defaults to 720 when left empty.

# Outputs
- `Pmu2mup`: probabilities for scattering in 3D from beam to beam. Matrix [n`_`direction x
    n`_`direction]
- `theta2beamW`: weight of each sub-beam within each beam. Matrix [n`_`beam x
    n`_`direction]
- `BeamWeight_relative`: relative weight of each sub-beam within each beam. It is the same
    as theta2beamW but normalized so that summing along the sub-beams gives 1 for each beam.
    Matrix [n`_`beam x n`_`direction]
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
                    verbose && print("Loading scattering-matrices from file: ", scattering_files[i1])
                    file = matopen(filename)
                    Pmu2mup = read(file, "Pmu2mup")
                    theta2beamW = read(file, "theta2beamW")
                    BeamWeight_relative = read(file, "BeamWeight_relative")
                    θ₁ = read(file, "theta1")
                    close(file)
                    found_them = 1
                    verbose && println(" ✅")
                    break
                end
            catch
            end
        end
    end
    if found_them == 0
        verbose && println("Could not find file with matching pitch-angle grid.")
        verbose && println("Starting to calculate the requested scattering-matrices.")

        Pmu2mup, theta2beamW, BeamWeight_relative, θ₁ = calculate_scattering_matrices(θ_lims, n_direction)

        # Save the results for future use
        filename = pkgdir(AURORA, "internal_data", "e_scattering",
                            string(length(θ_lims) - 1, "_streams_",
                            Dates.format(now(), "yyyymmdd-HHMMSS"),
                            ".mat"))
        file = matopen(filename, "w")
        write(file, "Pmu2mup", Pmu2mup)
        write(file, "theta2beamW", theta2beamW)
        write(file, "BeamWeight_relative", BeamWeight_relative)
        write(file, "theta1", θ₁)
        write(file, "theta_lims", Vector(θ_lims))
        write(file, "n_direction", n_direction)
        close(file)
    end

    return Pmu2mup, theta2beamW, BeamWeight_relative, θ₁
end






# Another 100x factor (0.2s vs 20s for the legacy version, still with θ_lims = 180:-10:0 and
# n_direction = 720
"""
    calculate_scattering_matrices(θ_lims, n_direction = 720)

Calculate the scattering matrices for given pitch-angle limits `θ_lims` of the electron
beams. Uses 720 directions by default for the start angle and the scattering angles,
equivalent to (1/4)° steps.

# Calling
`Pmu2mup, theta2beamW, BeamWeight_relative, θ₁ =  calculate_scattering_matrices(θ_lims, n_direction)`

# Inputs
- `θ_lims`: pitch-angle limits of the electron beams (e.g. 180:-10:0), where 180°
    corresponds to field aligned down, and 0° field aligned up.
- `n_direction`: number of directions or sub-beams to use for the discretized calculations
    of the scattering matrices. Defaults to 720 when left empty.

# Outputs
- `Pmu2mup`: probabilities for scattering in 3D from beam to beam. Matrix [n\\_direction x * n\\_direction x n\\_beam]
- `theta2beamW`: weight of each sub-beam within each beam. Matrix [n\\_beam x n\\_direction]
- `BeamWeight_relative`: relative weight of each sub-beam within each beam. It is the same
    as theta2beamW but normalized so that summing along the sub-beams gives 1 for each beam.
    Matrix [n\\_beam x n\\_direction]
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

    Pmu2mup = Array{Float64}(undef, length(θ₀), length(θ₁), length(θ_lims) - 1)
    theta2beamW = Array{Float64}(undef, length(θ_lims) - 1, length(θ₀))

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
                Pmu2mup[i0, i1, iμ] = 1 / π *
                                  (asin(clamp((z₂ - ct₀[i0]ct₁[i1]) / (st₀[i0]st₁[i1]), -1, 1)) -
                                   asin(clamp((z₁ - ct₀[i0]ct₁[i1]) / (st₀[i0]st₁[i1]), -1, 1)))

            end
        end
        next!(p)
    end

    # Normalize so that all sum(Pmu2mup[i, j, :], dims=3) = 1
    Pmu2mup = Pmu2mup ./ repeat(sum(Pmu2mup, dims = 3), outer = (1, 1, size(Pmu2mup, 3)))
    # Calculate the beam weight of each subdivision
    for iμ in (length(μ_lims) - 1):-1:1
        theta2beamW[iμ, :] .= abs.(sin.(θ₀)) .* (μ_lims[iμ] .< cos.(θ₀) .< μ_lims[iμ + 1])
    end
    # Normalize so that all sum(BeamWeight_relative, dims=2) = 1
    beam_sums = sum(theta2beamW, dims=2)
    if any(iszero, beam_sums)
        error(
            "Some pitch-angle beams have no sub-beam coverage. " *
            "This typically happens when θ_lims doesn't span the full range from 180° to 0°. " *
            "Make sure θ_lims includes both 180° and 0°."
        )
    end
    BeamWeight_relative = theta2beamW ./ repeat(beam_sums, 1, size(theta2beamW, 2))

    return Pmu2mup, theta2beamW, BeamWeight_relative, θ₁
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
