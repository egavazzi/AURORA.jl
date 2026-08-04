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

    # Production path: MSISDensity reads the file into a VectorDensity, which each species
    # samples on the grid (as in initialize!(model)).
    for species in (:N2, :O2, :O)
        source = MSISDensity(msis_file, species)
        @test source isa VectorDensity
        n = source(z)
        @test !any(isnan.(n))
        @test !any(isinf.(n))
        @test !any(n .< 0)
    end
end


@testitem "Electron densities" begin
    z = make_altitude_grid(50, 800)
    iri_file = find_iri_file(; verbose=false)
    # Production path: read_iri_file reads the file into an ElectronProfile, sampled on the grid
    source = read_iri_file(iri_file)
    @test source isa ElectronProfile
    ne, Te = source(z)

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
    ne, Te = read_iri_file(iri_file)(z)

    @test !any(isnan.(ne))
    @test !any(isinf.(ne))
    @test !any(ne .< 0)

    @test !any(isnan.(Te))
    @test !any(isinf.(Te))
    @test !any(Te .< 0)
end


@testitem "ElectronProfile sources (CCMC, legacy file, BYO)" begin
    z = make_altitude_grid(50, 800)

    # CCMC ModelWeb IRI export: variable preamble, Ne in cm⁻³, -1 sentinels
    ccmc_file = joinpath(@__DIR__, "test_data", "iri_ccmc_output.txt")
    p = read_ccmc_iri(ccmc_file)
    @test p isa ElectronProfile
    @test occursin("CCMC IRI", p.source)
    @test all(p.ne .> 0) && all(p.Te .> 0)        # sentinel rows dropped
    @test 1e11 < maximum(p.ne) < 1e13             # m⁻³ (cm⁻³ × 1e6), F-peak ~1e12
    ne, Te = p(z)
    @test all(isfinite, ne) && all(ne .> 0)
    @test all(isfinite, Te) && all(Te .> 0)

    # Legacy AURORA IRI file: read_iri_file matches the Ionosphere-via-path result exactly
    iri_file = joinpath(@__DIR__, "test_data", "iri_20051008-2200_70N-19E.txt")
    pl = read_iri_file(iri_file)
    @test occursin("IRI file", pl.source)
    ne1, Te1 = pl(z)
    iono = Ionosphere("label", iri_file, z)       # path → to_electron_source → read_iri_file
    @test iono.electron_source isa ElectronProfile
    @test ne1 ≈ iono.ne rtol=1e-9
    @test Te1 ≈ iono.Te rtol=1e-9

    # Bring-your-own vectors
    e = ElectronProfile([40e3, 100e3, 300e3, 850e3], [1e10, 1e11, 8e11, 5e10],
                        [250.0, 400.0, 1000.0, 1400.0]; source="byo")
    ne2, Te2 = e(z)
    @test length(ne2) == length(z) && all(ne2 .> 0) && all(Te2 .> 0)

    # A non-CCMC file is rejected with a clear error
    @test_throws ArgumentError read_ccmc_iri(iri_file)
end

@testitem "NeutralProfile sources (CCMC, legacy file)" begin
    z = make_altitude_grid(100, 600)

    # CCMC ModelWeb NRLMSIS export: preamble note, densities in cm⁻³, 9.999E-38 sentinels
    ccmc_file = joinpath(@__DIR__, "test_data", "msis_ccmc_output.txt")
    p = read_ccmc_msis(ccmc_file)
    @test p isa NeutralProfile
    @test occursin("CCMC NRLMSIS", p.source)
    @test all(haskey(p, s) for s in (:N2, :O2, :O))

    n_N2 = p[:N2](z)
    @test all(isfinite, n_N2) && all(n_N2 .> 0)
    @test issorted(n_N2; rev=true)                # N₂ falls off with altitude
    @test 1e14 < n_N2[1] < 1e19                   # m⁻³ (cm⁻³ × 1e6) near 100 km

    # O is sentinel-valued below ~100 km in this export; those levels are dropped, so the
    # O source starts higher than the N₂ one rather than carrying 9.999E-38 garbage
    @test minimum(p[:O].h) > minimum(p[:N2].h)
    @test all(p[:O].n .> 1e-37)

    # Legacy AURORA MSIS file: read_msis_file agrees with the per-species MSISDensity
    msis_file = find_msis_file(; verbose=false)
    np = read_msis_file(msis_file)
    @test occursin("MSIS file", np.source)
    for species in (:N2, :O2, :O)
        @test np[species](z) ≈ MSISDensity(msis_file, species)(z) rtol=1e-12
    end

    # A non-CCMC file is rejected with a clear error
    @test_throws ArgumentError read_ccmc_msis(msis_file)
end

@testitem "CCMC readers resolve columns by header name" begin
    mktempdir() do dir
        # An export whose columns sit in a different order than the reference file, and
        # which reports only some of the species, must still be read correctly.
        msis = joinpath(dir, "shuffled_msis.txt")
        write(msis, """
            Heit(km) N2den(cm-3) Oden(cm-3) O2den(cm-3)
            100.0    6.884E+12   4.793E+11  1.616E+12
            150.0    2.690E+10   7.917E+09  1.965E+09
            200.0    3.127E+09   2.100E+09  1.338E+08
            """)
        p = read_ccmc_msis(msis)
        @test sort(collect(keys(p))) == [:N2, :O, :O2]
        @test p[:N2].h ≈ [100e3, 150e3, 200e3]
        @test p[:N2].n ≈ [6.884e18, 2.690e16, 3.127e15]
        @test p[:O].n  ≈ [4.793e17, 7.917e15, 2.100e15]

        # Same for the IRI table
        iri = joinpath(dir, "shuffled_iri.txt")
        write(iri, """
              km   Te/K  Ne/cm-3
             100.0  159    48077
             200.0 1439   250319
             300.0 1740  1011143
            """)
        e = read_ccmc_iri(iri)
        @test e.h  ≈ [100e3, 200e3, 300e3]
        @test e.ne ≈ [4.8077e10, 2.50319e11, 1.011143e12]
        @test e.Te ≈ [159.0, 1439.0, 1740.0]

        # A required column that is missing is reported, with the header listed
        no_altitude = joinpath(dir, "no_altitude.txt")
        write(no_altitude, "N2den(cm-3) Oden(cm-3)\n6.884E+12 4.793E+11\n")
        @test_throws "Columns found: N2den(cm-3), Oden(cm-3)" read_ccmc_msis(no_altitude)

        no_Te = joinpath(dir, "no_Te.txt")
        write(no_Te, "km Ne/cm-3\n100.0 48077\n")
        @test_throws "no \"Te/K\" column" read_ccmc_iri(no_Te)

        # A header with no recognised species column is reported too
        no_species = joinpath(dir, "no_species.txt")
        write(no_species, "Heit(km) N2den(m-3)\n100.0 6.884E+18\n")
        @test_throws "no known species density column" read_ccmc_msis(no_species)
    end
end
