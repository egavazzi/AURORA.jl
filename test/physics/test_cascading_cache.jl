@testitem "CascadingCache provides spectra accessors" begin
    energy_grid = AURORA.EnergyGrid(100)
    cache = AURORA.CascadingCache()

    AURORA.load_or_compute_cascading!(cache, energy_grid)
    i_primary = searchsortedlast(cache[1].E_edges, 40.0)
    secondary = AURORA.secondary_spectrum(cache[1], i_primary, 15.581)
    primary = AURORA.primary_spectrum(cache[1], i_primary, 15.581)

    @test length(cache) == 3
    @test length(secondary) == energy_grid.n
    @test length(primary) == energy_grid.n
    @test sum(secondary) > 0
    @test sum(primary) >= 0
    @test !isempty(cache[1].E_edges)
    @test !isempty(cache[1].ionization_thresholds)
    @test size(cache[1].secondary_transfer_matrix, 1) == energy_grid.n
    @test size(cache[1].secondary_transfer_matrix, 2) == energy_grid.n
end
