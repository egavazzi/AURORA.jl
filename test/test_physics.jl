@testitem "ScatteringData construction" begin
    θ_lims = 180:-30:0
    sd = AURORA.ScatteringData(θ_lims)

    @test sd isa AURORA.ScatteringData
    n_beams = length(θ_lims) - 1
    @test size(sd.Ω_subbeam_relative, 1) == n_beams
    @test size(sd.P_scatter, 1) == size(sd.P_scatter, 2)
    @test length(sd.Ω_beam) == n_beams
    @test length(sd.θ_scatter) > 0
end

@testitem "CrossSectionData construction" begin
    energy_grid = AURORA.EnergyGrid(500)
    cs = AURORA.CrossSectionData(energy_grid)

    @test cs isa AURORA.CrossSectionData
    @test size(cs.σ_neutrals.σ_N2, 2) == energy_grid.n
    @test size(cs.σ_neutrals.σ_O2, 2) == energy_grid.n
    @test size(cs.σ_neutrals.σ_O, 2) == energy_grid.n
    @test size(cs.collision_levels.N2_levels, 2) == 2
    @test size(cs.collision_levels.O2_levels, 2) == 2
    @test size(cs.collision_levels.O_levels, 2) == 2
end
