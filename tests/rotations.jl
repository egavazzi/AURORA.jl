using Rotations
using GLMakie
using AURORA
using ProgressMeter

θ_lims = 180:-10:0
μ_lims = cosd.(θ_lims)

n_dirs = 720

function rotating(θ_lims, n_dirs)
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

    return Pmu2mup, B, theta2beamW
end

##
# using BenchmarkTools
# @profview Pmu2mup = rotating(θ_lims, 180)
# @profview_allocs Pmu2mup = rotating(θ_lims, n_dirs)
# @btime Pmu2mup = rotating(θ_lims, n_dirs);

Pmu2mup, B, theta2beamW = rotating(θ_lims, n_dirs);

# all(sum(B, dims=3) .== n_dirs)

θ₀ = acos.(mu_avg(0:(180/n_dirs):180))
# a = vec(sum(sin.(θ₀) .* Pmu2mup, dims=(1, 2))); a = a ./ sum(a) * 2pi

a = similar(Pmu2mup);
for i in axes(theta2beamW, 1)
    a[:, :, i] = theta2beamW[i, :] .* Pmu2mup[:, :, i]
end
a = vec(sum(a, dims=(1, 2))); a = a ./ sum(a) * 4pi
(a .- b)


##
test = sum(theta2beamW, dims=2); test = test / sum(test) * 4pi
b = beam_weight(θ_lims)

(test .- b) ./ b

##
μ_lims = cosd.(θ_lims)

bla = reduce(vcat, bla')

for iμ in (length(μ_lims) - 1):-1:1
    cache1 = (μ_lims[iμ] .< bla[:, 3] .< μ_lims[iμ + 1])
    println(sum(cache1))
end
# for iμ in (length(μ_lims) - 1):-1:1
#     cache1 = (μ_lims[iμ] .< bla[(end-180:end), 3] .< μ_lims[iμ + 1])
#     println(sum(cache1))
# end



(μ_lims[1] .< bla[179*180 .+ (1:181), 3] .< μ_lims[1 + 1])
##

for i in 1:length(θ_lims) - 1
    println(sum(Pmu2mup[:, :, i]) / sum(Pmu2mup))
end

for i in 1:length(θ_lims) - 1
    println(sum(Pmu2mup[:, :, i]))
end

# for i in 1:Int((length(θ_lims) - 1)/2)
#     println(sum(Pmu2mup[:, :, i] - reverse(Pmu2mup[:, :, i], dims=1)))
# end

##
using Test
@testset begin
    for i in axes(Pmu2mup, 1)
        @test sum(Pmu2mup[i, :, :]) ≈ size(Pmu2mup, 2)
    end
end
@testset begin
    for i in axes(Pmu2mup, 2)
        @test sum(Pmu2mup[:, i, :]) ≈ size(Pmu2mup, 1)
    end
end
@testset begin
    for i in axes(Pmu2mup, 1), j in axes(Pmu2mup, 2)
        @test sum(Pmu2mup[i, j, :]) ≈ 1
    end
end
##


# a = vec(sum(Pmu2mup, dims=(1, 2)) ./ sum(Pmu2mup))
b = beam_weight(θ_lims);

θ₀ = acos.(mu_avg(0:(180/n_dirs):180))
a = vec(sum(sin.(θ₀) .* Pmu2mup, dims=(1, 2)) ./ sum(Pmu2mup)); a = a / sum(a) * 4pi

((a .- b) ./ b) # relative diff with beam weight
a .- b




θ₀ = acos.(mu_avg(0:(180/n_dirs):180))
θ₁ = deg2rad.(0:(180/n_dirs):360); θ₁ = θ₁[1:end-1]
ϕ = deg2rad.(0:(180/n_dirs):360); #ϕ = ϕ[1:end-1]
0.00025365486373079074
-0.00025365486373085666
-0.00025365486373080115
0.0002536548637308185


θ₀ = acos.(mu_avg(0:(180/n_dirs):180))
θ₁ = deg2rad.(0:(180/n_dirs):360); θ₁ = θ₁[1:end-1]
ϕ = deg2rad.(0:(180/n_dirs):360); ϕ = ϕ[1:end-1]
2.947773417430813e-5
-2.9477734174387926e-5
-2.9477734174332415e-5
2.9477734174335884e-5





θ₀ = acos.(mu_avg(0:(180/n_dirs):180))
θ₁ = deg2rad.(0:(180/n_dirs):360); θ₁ = θ₁[1:end-1]
ϕ = deg2rad.(0:(180/n_dirs):180); ϕ = ϕ[1:end-1]
2.947773417430813e-5
-2.9477734174276904e-5
-2.9477734174221393e-5
2.9477734174335884e-5


θ₀ = acos.(mu_avg(0:(180/n_dirs):180))
θ₁ = deg2rad.(0:(180/n_dirs):360); θ₁ = θ₁[1:end-1]
ϕ = deg2rad.(0:(180/n_dirs):180); #ϕ = ϕ[1:end-1]
0.0004765934456101642
-0.00047659344561024053
-0.000476593445610185
0.00047659344561019196






θ₀ = acos.(mu_avg(0:(180/n_dirs):180))
θ₁ = deg2rad.(0:(180/n_dirs):180); θ₁ = θ₁[1:end-1]
ϕ = deg2rad.(0:(180/n_dirs):360); ϕ = ϕ[1:end-1]
2.9618560976388958e-5
-2.961856097638549e-5
-2.9618560976329977e-5
2.9618560976416713e-5


θ₀ = acos.(mu_avg(0:(180/n_dirs):180))
θ₁ = deg2rad.(0:(180/n_dirs):180); #θ₁ = θ₁[1:end-1]
ϕ = deg2rad.(0:(180/n_dirs):360); ϕ = ϕ[1:end-1]
2.975860972984079e-5
-3.3031944386752166e-5
-2.6485275073040437e-5
2.9758609729868546e-5

θ₀ = acos.(mu_avg(0:(180/n_dirs):180))
θ₁ = deg2rad.(0:(180/n_dirs):360); #θ₁ = θ₁[1:end-1]
ϕ = deg2rad.(0:(180/n_dirs):360); ϕ = ϕ[1:end-1]
2.9618560976388958e-5
-2.961856097638549e-5
-2.9618560976329977e-5
2.9618560976416713e-5


θ₀ = acos.(mu_avg(0:(180/n_dirs):180))
θ₁ = deg2rad.(0:(180/n_dirs):180); θ₁ = θ₁[1:end-1]
ϕ = deg2rad.(0:(180/n_dirs):180); ϕ = ϕ[1:end-1]
0.0006342050540026171
0.0011380875791615064
-0.0011970430475102267
-0.0005752495856539835
