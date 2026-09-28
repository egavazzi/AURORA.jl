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

@testitem "Ionizing levels match the default cascading specs" begin
    specs = Dict("N2" => AURORA.DefaultCascadingSpecN2(),
                 "O2" => AURORA.DefaultCascadingSpecO2(),
                 "O"  => AURORA.DefaultCascadingSpecO())

    for (species, spec) in specs
        E_levels = AURORA.load_excitation_threshold_for(species)
        ionizing = findall(>(0), E_levels[:, 2])
        @test E_levels[ionizing, 1] == spec.ionization_thresholds
        @test E_levels[ionizing, 2] == spec.n_secondaries
    end

    # O⁺(2s2p⁴ ⁴P) is a single ionization with one secondary electron.
    names = AURORA.get_level_names("O")
    O_levels = AURORA.load_excitation_threshold_for("O")
    i_ionion = only(findall(==("ionion"), names))
    @test O_levels[i_ionion, 2] == 1
end
