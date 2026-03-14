@testitem "Ionosphere construction" begin
    h_atm = make_altitude_grid(50, 800)
    msis_file = find_msis_file()
    iri_file = find_iri_file()
    iono = Ionosphere(msis_file, iri_file, h_atm)

    @test iono isa Ionosphere
    @test length(iono.nN2) == length(h_atm)
    @test length(iono.nO2) == length(h_atm)
    @test length(iono.nO) == length(h_atm)
    @test length(iono.Tn) == length(h_atm)
    @test length(iono.Te) == length(h_atm)
    @test length(iono.ne) == length(h_atm)

    for field in (iono.nN2, iono.nO2, iono.nO, iono.Tn, iono.Te, iono.ne)
        @test !any(isnan.(field))
        @test !any(isinf.(field))
        @test !any(field .< 0)
    end

    nn = AURORA.n_neutrals(iono)
    @test nn.nN2 === iono.nN2
    @test nn.nO2 === iono.nO2
    @test nn.nO === iono.nO
end
