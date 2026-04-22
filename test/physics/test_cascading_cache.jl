@testitem "CascadingCache provides spectra accessors" begin
    energy_grid = AURORA.EnergyGrid(100)
    E_grid = energy_grid.E_edges[1:end-1]
    dE = energy_grid.ΔE
    cache = AURORA.CascadingCache()

    secondary = AURORA.secondary_spectrum(cache[1], E_grid, dE, 40.0, 15.581)
    primary = AURORA.primary_spectrum(cache[1], E_grid, dE, 40.0, 15.581)

    @test length(cache) == 3
    @test length(secondary) == length(E_grid)
    @test length(primary) == length(E_grid)
    @test sum(secondary) > 0
    @test sum(primary) >= 0
    @test !isempty(cache[1].E_grid_for_Q)
    @test !isempty(cache[1].ionization_thresholds)
end
