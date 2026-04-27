@testitem "CascadingCache provides spectra accessors" begin
    energy_grid = AURORA.EnergyGrid(100)
    cache = AURORA.CascadingCache()

    AURORA.load_or_compute_cascading!(cache, energy_grid)
    secondary = AURORA.secondary_spectrum(cache[1], energy_grid.E_centers, 40.0, 15.581)
    primary = AURORA.primary_spectrum(cache[1], 40.0, 15.581)

    @test length(cache) == 3
    @test length(secondary) == energy_grid.n
    @test length(primary) == energy_grid.n
    @test sum(secondary) > 0
    @test sum(primary) >= 0
    @test !isempty(cache[1].E_edges_for_Q)
    @test !isempty(cache[1].ionization_thresholds)
end
