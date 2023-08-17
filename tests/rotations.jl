using Rotations
using GLMakie
using AURORA
using ProgressMeter

θ_lims = 180:-10:0
n_dirs = 720

function rotating(θ_lims, n_dirs=720)
    μ_lims = cosd.(θ_lims)

    θ₀ = acos.(mu_avg(0:(180/n_dirs):180))
    θ₁ = deg2rad.(0:(180/n_dirs):360); θ₁ = θ₁[1:end-1] # θ₁ = θ₁[[1:n_dirs-1; n_dirs+1:end]]
    ϕ = deg2rad.(0:(180/n_dirs):180); ϕ = ϕ[1:end-1]

    μ = Vector{Float64}(undef, length(ϕ))
    B = Array{Float64}(undef, length(θ₀), length(θ₁), length(θ_lims) - 1)

    cache = Vector{Float64}(undef, length(ϕ))

    p = Progress(length(θ₀); desc=string("Calculating scattering matrices "), color=:blue);

    # bla = Vector[]
    for i0 in length(θ₀):-1:1
        e0 = [sin(θ₀[i0]), 0, cos(θ₀[i0])]

        for i1 in length(θ₁):-1:1
            rotation_axis = [0, 1, 0]
            rotation_angle = θ₁[i1]
            rotation_matrix = AngleAxis(rotation_angle, rotation_axis...)
            e1 = rotation_matrix * e0

            # push!(bla, e1)
            for iϕ in length(ϕ):-1:1
                rotation_axis = e0
                rotation_angle = ϕ[iϕ]
                rotation_matrix = AngleAxis(rotation_angle, rotation_axis...)
                es = rotation_matrix * e1
                μ[iϕ] = es[3]
            end

            for iμ in (length(μ_lims) - 1):-1:1
                cache .= (μ_lims[iμ] .< μ .< μ_lims[iμ + 1])
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

    return Pmu2mup, theta2beamW
end
