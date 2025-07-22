# Define some parameters for the tests
θ_lims = 180:-10:0
n_dirs = 720
Pmu2mup, theta2beamW, BeamWeight_relative, θ₁ = AURORA.calculate_scattering_matrices(θ_lims, n_dirs);
BeamW = beam_weight(θ_lims);

# check if θ_lims is symmetric
check = false
if iseven(length(θ_lims) - 1)
    n_θhalf = Int((length(θ_lims) - 1) / 2)
    check = all(θ_lims[1:n_θhalf + 1] .+ θ_lims[end:-1:n_θhalf + 1] .== 180)
end


@testset "Rotation matrices" begin
    @testset "BeamWeight_relative" begin
        # check normalization
        @test all(sum(BeamWeight_relative, dims=2) .≈ 1)
    end
    @testset "theta2beamW" begin
        # summing and normalizing theta2beamW along the beams should give the beam weights
        BW = sum(theta2beamW, dims=2);
        BW = BW / sum(BW) * 4pi;
        @test BW ≈ BeamW
    end

    @testset "Pmu2mup" begin
        # check normalization (for any i and j, sum(Pmu2mup[i, j, :], dims=3) = 1)
        @test all(sum(Pmu2mup, dims = 3) .≈ 1)
        if check
            # if θ_lims is symmetric, check that sum(Pmu2mup, dims=(1, 2)) is symmetric
            for i in 1:n_θhalf
                @test isapprox(sum(Pmu2mup[:, :, i]), sum(Pmu2mup[:, :, end + 1 - i]))
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

        B2B_fragment = AURORA.prepare_beams2beams(BeamWeight_relative, Pmu2mup);

        iE = 400 # ~3730 eV (just to pick one)
        for i in 1:3
            phase_fcn_e = convert_phase_fcn_to_3D(phase_fcn[i][1][:, iE], θ₁);
            phase_fcn_i = convert_phase_fcn_to_3D(phase_fcn[i][2][:, iE], θ₁);
            B2B_elastic = beams2beams(phase_fcn_e, B2B_fragment);
            B2B_inelastic = beams2beams(phase_fcn_i, B2B_fragment);

            @test all(sum(B2B_elastic, dims = 1) .≈ 1)
            @test all(sum(B2B_inelastic, dims = 1) .≈ 1)
        end
    end
end
