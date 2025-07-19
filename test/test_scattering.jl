using AURORA
using Test

θ_lims = 180:-10:0
n_dirs = 720
Pmu2mup, theta2beamW, BeamWeight_relative, θ₁ = rotating(θ_lims, n_dirs);
BeamW = beam_weight(θ_lims);

# check if θ_lims is symmetric
check = false
if iseven(length(θ_lims) - 1)
    n_θhalf = Int((length(θ_lims) - 1) / 2)
    check = all(θ_lims[1:n_θhalf + 1] .+ θ_lims[end:-1:n_θhalf + 1] .== 180)
end


@testset "Rotation matrices" begin
    # summing and normalizing theta2beamW along the beams should be equal to the beam weights
    @testset "theta2beamW" begin
        BW = sum(theta2beamW, dims=2);
        BW = BW / sum(BW) * 4pi;
        @test BW ≈ BeamW
    end

    # if θ_lims is symmetric, sum(Pmu2mup, dims=(1, 2)) should be symmetric
    if check
        @testset "Pmu2mup" begin
            for i in 1:n_θhalf
            @test (sum(Pmu2mup[:, :, i]) - sum(Pmu2mup[:, :, end + 1 - i])) ≈ 0
            end
        end
    end
end


# if θ_lims is symmetric, the B2B matrices should be symmetric
if check
    @testset "B2B matrices" begin
        E, _ = AURORA.make_energy_grid(7000)
        phaseN2e, phaseN2i = phase_fcn_N2(θ₁, E);
        phaseO2e, phaseO2i = phase_fcn_O2(θ₁, E);
        phaseOe, phaseOi = phase_fcn_O(θ₁, E);
        phase_fcn = ((phaseN2e, phaseN2i), (phaseO2e, phaseO2i), (phaseOe, phaseOi));

        iE = 400
        for i in 1:3
            phase_fcn_e = convert_phase_fcn_to_3D(phase_fcn[i][1][:, iE], θ₁);
            phase_fcn_i = convert_phase_fcn_to_3D(phase_fcn[i][2][:, iE], θ₁);
            B2B_elastic = beams2beams(phase_fcn_e, Pmu2mup, BeamWeight_relative);
            B2B_inelastic = beams2beams(phase_fcn_i, Pmu2mup, BeamWeight_relative);

            @test all(sum(B2B_elastic, dims = 1) .≈ 1)
            @test all(sum(B2B_inelastic, dims = 1) .≈ 1)
        end
    end
end
