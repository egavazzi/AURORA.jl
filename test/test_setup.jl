@testitem "Neutral densities" begin
    h_atm = make_altitude_grid(50, 800)
    msis_file = find_msis_file()

    @test_nowarn n_neutrals, Tn = AURORA.load_neutral_densities(msis_file, h_atm)

    n_neutrals, Tn = AURORA.load_neutral_densities(msis_file, h_atm)
    for i in eachindex(n_neutrals)
        @test !any(isnan.(n_neutrals[i]))
        @test !any(isinf.(n_neutrals[i]))
        @test !any(n_neutrals[i] .< 0)
    end

    @test !any(isnan.(Tn))
    @test !any(isinf.(Tn))
    @test !any(Tn .< 0)
end


@testitem "Electron densities" begin
    h_atm = make_altitude_grid(50, 800)
    iri_file = find_iri_file()
    ne, Te = AURORA.load_electron_densities(iri_file, h_atm)

    @test !any(isnan.(ne))
    @test !any(isinf.(ne))
    @test !any(ne .< 0)

    @test !any(isnan.(Te))
    @test !any(isinf.(Te))
    @test !any(Te .< 0)
end
