@testitem "Cross-section loading functions" begin
    energy_grid = AURORA.EnergyGrid(500)

    σ_N2 = AURORA.get_cross_section("N2", energy_grid.E_centers)
    σ_O2 = AURORA.get_cross_section("O2", energy_grid.E_centers)
    σ_O  = AURORA.get_cross_section("O",  energy_grid.E_centers)

    @test size(σ_N2, 2) == energy_grid.n
    @test size(σ_O2, 2) == energy_grid.n
    @test size(σ_O,  2) == energy_grid.n

    N2_levels = AURORA.load_excitation_threshold_for("N2")
    O2_levels = AURORA.load_excitation_threshold_for("O2")
    O_levels  = AURORA.load_excitation_threshold_for("O")

    @test size(N2_levels, 2) == 2
    @test size(O2_levels, 2) == 2
    @test size(O_levels,  2) == 2
end
