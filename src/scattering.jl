using Rotations
using AURORA
using ProgressMeter
using MAT

# 144x faster than the Matlab code (24s instead of 1h for θ_lims = 180:-10:0 and
# n_direction = 720)
"""
    rotating(θ_lims, n_direction=720)

Calculate the scattering matrices for given pitch-angle limits `θ_lims` of the electron
beams.

# Calling
`Pmu2mup, theta2beamW, BeamWeight_relative, θ₁ =  rotating(θ_lims, n_direction)`

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
function rotating(θ_lims, n_direction=720)
    μ_lims = cosd.(θ_lims)

    θ₀ = Vector(0:(180/n_direction):180); θ₀ = deg2rad.((θ₀[1:end - 1] .+ θ₀[2:end]) / 2)
    θ₁ = deg2rad.(0:(180/n_direction):360); θ₁ = θ₁[1:end-1]
    ϕ = deg2rad.(0:(180/n_direction):360); ϕ = ϕ[1:end-1]

    μ = Vector{Float64}(undef, length(ϕ))
    B = Array{Float64}(undef, length(θ₀), length(θ₁), length(θ_lims) - 1)

    cache = Vector{Float64}(undef, length(ϕ))

    p = Progress(length(θ₀); desc=string("Calculating scattering matrices "), color=:blue);

    rotation_axis1 = [0, 1, 0]
    # e0 = zeros(3)
    for i0 in length(θ₀):-1:1
        e0 = [sin(θ₀[i0]), 0, cos(θ₀[i0])]
        # e0[1] = sin(θ₀[i0])
        # e0[2] = 0
        # e0[3] = cos(θ₀[i0])

        for i1 in length(θ₁):-1:1
            rotation_angle = θ₁[i1]
            rotation_matrix = AngleAxis(rotation_angle, rotation_axis1[1], rotation_axis1[2], rotation_axis1[3])
            e1 = rotation_matrix * e0

            for iϕ in length(ϕ):-1:1
                rotation_axis = e0
                rotation_angle = ϕ[iϕ]
                rotation_matrix = AngleAxis(rotation_angle, rotation_axis[1], rotation_axis[2], rotation_axis[3])
                es = rotation_matrix * e1
                μ[iϕ] = es[3]
            end

            for iμ in (length(μ_lims) - 1):-1:1
                cache .= μ_lims[iμ] .< μ .< μ_lims[iμ + 1]
                B[i0, i1, iμ] = sum(cache)
            end
        end
        next!(p)
    end

    Pmu2mup = B ./ repeat(sum(B, dims=3), outer = (1, 1, size(B, 3)))
    # Pmu2mup[1, 1, :] .= Pmu2mup[2, 1, :] ./ 2 .+ Pmu2mup[1, 2, :] ./ 2
    # Pmu2mup[end, end, :] .= Pmu2mup[end - 1, end, :] ./ 2 .+ Pmu2mup[end, end - 1, :] ./ 2

    theta2beamW = Array{Float64}(undef, length(θ_lims) - 1, length(θ₀))
    for iμ in (length(μ_lims) - 1):-1:1
        theta2beamW[iμ, :] .= abs.(sin.(θ₀)) .* (μ_lims[iμ] .< cos.(θ₀) .< μ_lims[iμ + 1])
    end
    BeamWeight_relative = theta2beamW ./ repeat(sum(theta2beamW, dims=2), 1, size(theta2beamW, 2))

    return Pmu2mup, theta2beamW, BeamWeight_relative, θ₁
end


"""
    load_scattering_matrices(θ_lims, n_direction=720)

Look for scattering matrices that match the pitch-angle limits `θ_lims` and the number
of direction/sub-beams `n_direction`. If a file is found, the scattering matrices are
directly loaded. Otherwise, they are calculated and saved to a file.

# Calling
`Pmu2mup, theta2beamW, BeamWeight_relative, θ₁ =  load_scattering_matrices(θ_lims, n_direction)`

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
function load_scattering_matrices(θ_lims, n_direction=720)
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
                    print("Loading scattering-matrices from file: ", scattering_files[i1])
                    file = matopen(filename)
                    Pmu2mup = read(file, "Pmu2mup")
                    theta2beamW = read(file, "theta2beamW")
                    BeamWeight_relative = read(file, "BeamWeight_relative")
                    θ₁ = read(file, "theta1")
                    close(file)
                    found_them = 1
                    println(" ✅")
                    break
                end
            catch
            end
        end
    end
    if found_them == 0
        println("Could not find file with matching pitch-angle grid.")
        println("Starting to calculate the requested scattering-matrices.")

        Pmu2mup, theta2beamW, BeamWeight_relative, θ₁ = rotating(θ_lims, n_direction)

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
