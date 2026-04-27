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
