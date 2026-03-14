@testitem "ScatteringData construction" begin
    θ_lims = 180:-30:0
    sd = ScatteringData(θ_lims)

    @test sd isa ScatteringData
    n_beams = length(θ_lims) - 1
    @test size(sd.BeamWeight_relative, 1) == n_beams
    @test size(sd.Pmu2mup, 1) == size(sd.Pmu2mup, 2)
    @test length(sd.BeamWeight) == n_beams
    @test length(sd.theta1) > 0
end

@testitem "CrossSectionData construction" begin
    E, dE = make_energy_grid(500)
    cs = CrossSectionData(E, dE)

    @test cs isa CrossSectionData
    @test size(cs.σ_neutrals.σ_N2, 2) == length(E)
    @test size(cs.σ_neutrals.σ_O2, 2) == length(E)
    @test size(cs.σ_neutrals.σ_O, 2) == length(E)
    @test size(cs.E_levels_neutrals.N2_levels, 2) == 2
    @test size(cs.E_levels_neutrals.O2_levels, 2) == 2
    @test size(cs.E_levels_neutrals.O_levels, 2) == 2
end
