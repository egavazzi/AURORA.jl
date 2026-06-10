@testitem "Ionosphere construction" begin
    z = make_altitude_grid(50, 800)
    msis_file = find_msis_file(; verbose=false)
    iri_file = find_iri_file(; verbose=false)
    iono = Ionosphere(msis_file, iri_file, z)

    @test iono isa Ionosphere
    @test length(iono.Te) == length(z)
    @test length(iono.ne) == length(z)

    for field in (iono.Te, iono.ne)
        @test !any(isnan.(field))
        @test !any(isinf.(field))
        @test !any(field .< 0)
    end
end


@testitem "Neutral densities" begin
    z = make_altitude_grid(50, 800)
    msis_file = find_msis_file(; verbose=false)

    # Production path: each species samples its MSIS density profile on the grid, then the
    # top boundary is tapered to zero by apply_density_boundary! (as in initialize!(model)).
    for species in (:N2, :O2, :O)
        n = AURORA.load_msis_density(msis_file, species, z)
        @test_nowarn AURORA.apply_density_boundary!(n)
        @test !any(isnan.(n))
        @test !any(isinf.(n))
        @test !any(n .< 0)
        @test all(iszero, n[end-2:end])   # tapered to zero at the top
    end
end


@testitem "Electron densities" begin
    z = make_altitude_grid(50, 800)
    iri_file = find_iri_file(; verbose=false)
    ne, Te = AURORA.load_electron_densities(iri_file, z)

    @test !any(isnan.(ne))
    @test !any(isinf.(ne))
    @test !any(ne .< 0)

    @test !any(isnan.(Te))
    @test !any(isinf.(Te))
    @test !any(Te .< 0)
end


@testitem "Recent electron densities from IRI" begin
    using Dates
    # Use a rolling date (1 year ago) to test the IRI coefficient files coverage
    past = today() - Year(1)
    z = make_altitude_grid(50, 800)
    iri_file = find_iri_file(; year = year(past), month = month(past), day = day(past),
                              hour = 12, minute = 0, verbose = false)
    ne, Te = AURORA.load_electron_densities(iri_file, z)

    @test !any(isnan.(ne))
    @test !any(isinf.(ne))
    @test !any(ne .< 0)

    @test !any(isnan.(Te))
    @test !any(isinf.(Te))
    @test !any(Te .< 0)
end
