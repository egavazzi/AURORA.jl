@testitem "Ionosphere construction" begin
    z = make_altitude_grid(50, 800)
    iri_file = find_iri_file(; verbose=false)
    iono = Ionosphere(iri_file, z)

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

    # Production path: read_msis_file reads the file into a DensityProfile per species, which
    # each species samples on the grid (as in initialize!(model)).
    neutrals = read_msis_file(msis_file)
    for species in (:N2, :O2, :O)
        source = neutrals[species]
        @test source isa DensityProfile
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
    @test occursin("CCMC IRI", p.origin)
    @test all(p.ne .> 0) && all(p.Te .> 0)        # sentinel rows dropped
    @test 1e11 < maximum(p.ne) < 1e13             # m⁻³ (cm⁻³ × 1e6), F-peak ~1e12
    ne, Te = p(z)
    @test all(isfinite, ne) && all(ne .> 0)
    @test all(isfinite, Te) && all(Te .> 0)

    # Legacy AURORA IRI file: read_iri_file matches the Ionosphere-via-path result exactly
    iri_file = joinpath(@__DIR__, "test_data", "iri_20051008-2200_70N-19E.txt")
    pl = read_iri_file(iri_file)
    @test occursin("IRI file", pl.origin)
    ne1, Te1 = pl(z)
    iono = Ionosphere(iri_file, z)
    @test iono.electron_source isa ElectronProfile
    @test ne1 ≈ iono.ne rtol=1e-9
    @test Te1 ≈ iono.Te rtol=1e-9

    # Bring-your-own vectors
    e = ElectronProfile([40e3, 100e3, 300e3, 850e3], [1e10, 1e11, 8e11, 5e10],
                        [250.0, 400.0, 1000.0, 1400.0]; origin="byo")
    ne2, Te2 = e(z)
    @test length(ne2) == length(z) && all(ne2 .> 0) && all(Te2 .> 0)

    # A non-CCMC file is rejected with a clear error
    @test_throws ArgumentError read_ccmc_iri(iri_file)
end

@testitem "NeutralAtmosphere sources (CCMC, legacy file)" begin
    z = make_altitude_grid(100, 600)

    # CCMC ModelWeb NRLMSIS export: preamble note, densities in cm⁻³, 9.999E-38 sentinels
    ccmc_file = joinpath(@__DIR__, "test_data", "msis_ccmc_output.txt")
    p = read_ccmc_msis(ccmc_file)
    @test p isa NeutralAtmosphere
    @test occursin("CCMC NRLMSIS", p.origin)
    @test all(haskey(p, s) for s in (:N2, :O2, :O))

    n_N2 = p[:N2](z)
    @test all(isfinite, n_N2) && all(n_N2 .> 0)
    @test issorted(n_N2; rev=true)                # N₂ falls off with altitude
    @test 1e14 < n_N2[1] < 1e19                   # m⁻³ (cm⁻³ × 1e6) near 100 km

    # O is sentinel-valued below ~100 km in this export; those levels are dropped, so the
    # O source starts higher than the N₂ one rather than carrying 9.999E-38 garbage
    @test minimum(p[:O].h) > minimum(p[:N2].h)
    @test all(p[:O].n .> 1e-37)

    # Legacy AURORA MSIS file: each species keeps the file's own native grid and values
    msis_file = find_msis_file(; verbose=false)
    np  = read_msis_file(msis_file)
    raw = AURORA.load_msis(msis_file)
    @test occursin("MSIS file", np.origin)
    for species in (:N2, :O2, :O)
        valid = AURORA.usable_levels(getproperty(raw.data, species))
        @test np[species].h ≈ raw.data.height_km[valid] .* 1e3 rtol=1e-12
        @test np[species].n ≈ getproperty(raw.data, species)[valid] rtol=1e-12
    end

    # A non-CCMC file is rejected with a clear error
    @test_throws ArgumentError read_ccmc_msis(msis_file)
end


@testitem "NeutralAtmosphere reports species it dropped" begin
    # A species the source carries but never usably reports is dropped at read time. Asking
    # for it later must say why, rather than looking like the species was never there.
    densities, dropped = AURORA.species_densities(
        [100e3, 200e3, 300e3],
        Dict(:N2 => [1e18, 1e17, 1e16], :NO => [NaN, NaN, NaN]),
        "test source")
    np = NeutralAtmosphere(densities; origin="test source", dropped)

    @test haskey(np, :N2)
    @test !haskey(np, :NO)
    @test np.dropped == [:NO]

    unreported = sprint(showerror, try np[:NO] catch e; e end)
    @test occursin("fewer than 2 usable levels", unreported)

    # A species that was never in the source at all gets no such hint
    absent = sprint(showerror, try np[:Ar] catch e; e end)
    @test !occursin("fewer than 2 usable levels", absent)
end

@testitem "run_msis returns a NeutralAtmosphere matching the file path" begin
    z = make_altitude_grid(100, 500)
    conditions = (; year = 2005, month = 10, day = 8, hour = 22, minute = 0,
                    lat = 69.58, lon = 19.23, height = 85:5:600)

    np = run_msis(; conditions..., verbose = false)
    @test np isa NeutralAtmosphere
    @test all(haskey(np, s) for s in (:N2, :O2, :O))
    @test occursin("NRLMSIS 2.1", np.origin)

    # MSIS reports no N below ~95 km (NaN); those levels are dropped for that species only,
    # so :N starts higher than :N2 instead of carrying NaN into the interpolation.
    @test minimum(np[:N].h) > minimum(np[:N2].h)
    for s in keys(np)
        @test all(isfinite, np[s].n) && all(np[s].n .> 0)
        @test all(isfinite, np[s](z))
    end

    # Same densities as going through a file. find_msis_file may return a file already in
    # the package store, and older stored files keep only ~7 significant digits.
    msis_file = find_msis_file(; conditions..., verbose = false)
    fp = read_msis_file(msis_file)
    for s in (:N2, :O2, :O)
        @test np[s](z) ≈ fp[s](z) rtol=1e-6
    end
end

@testitem "run_msis save_to writes a file that reads back identically" begin
    z = make_altitude_grid(100, 400)
    dir = mktempdir()
    np = run_msis(; year = 2005, month = 10, day = 8, hour = 22, minute = 0,
                    lat = 69.58, lon = 19.23, height = 85:5:600,
                    save_to = dir, verbose = false)

    files = readdir(dir)
    @test length(files) == 1
    @test startswith(files[1], "msis_") && endswith(files[1], ".txt")

    fp = read_msis_file(joinpath(dir, files[1]))
    @test fp isa NeutralAtmosphere
    @test sort(collect(keys(fp))) == sort(collect(keys(np)))
    for s in (:N2, :O2, :O)
        @test fp[s](z) ≈ np[s](z) rtol=1e-12
    end
end

@testitem "run_iri save_to writes a file that reads back identically" begin
    z = make_altitude_grid(100, 400)
    dir = mktempdir()
    p = run_iri(; year = 2005, month = 10, day = 8, hour = 22, minute = 0,
                  lat = 69.58, lon = 19.23, height = 85:5:600,
                  save_to = dir, verbose = false)

    files = readdir(dir)
    @test length(files) == 1
    @test startswith(files[1], "iri_") && endswith(files[1], ".txt")

    fp = read_iri_file(joinpath(dir, files[1]))
    @test fp isa ElectronProfile
    ne, Te = p(z)
    ne_file, Te_file = fp(z)
    @test ne_file ≈ ne rtol=1e-12
    @test Te_file ≈ Te rtol=1e-12
end

@testitem "Profile sources validate their vectors at construction" begin
    h = [100e3, 200e3, 300e3]

    # Length mismatch, too few levels, non-increasing altitudes, non-positive values
    @test_throws "must have the same length" DensityProfile(h, [1e18, 1e17])
    @test_throws "at least 2 altitude levels" DensityProfile([100e3], [1e18])
    @test_throws "strictly increasing" DensityProfile([300e3, 100e3, 200e3], [1e18, 1e17, 1e16])
    @test_throws "finite and strictly positive" DensityProfile(h, [1e18, 0.0, 1e16])
    @test_throws "finite and strictly positive" DensityProfile(h, [1e18, NaN, 1e16])

    @test_throws "must have the same length" ElectronProfile(h, [1e11, 1e12], [200.0, 300.0])
    @test_throws "finite and strictly positive" ElectronProfile(h, [1e11, -1.0, 1e12],
                                                                [200.0, 300.0, 400.0])
    @test_throws "finite and strictly positive" ElectronProfile(h, [1e11, 1e12, 1e12],
                                                                [200.0, 0.0, 400.0])

    # A valid profile still builds, and reports where it came from
    @test DensityProfile(h, [1e18, 1e17, 1e16]; origin = "ok") isa DensityProfile
    @test occursin("ok", sprint(show, DensityProfile(h, [1e18, 1e17, 1e16]; origin = "ok")))
end

@testitem "Sampling outside a profile's native range warns" begin
    h = [100e3, 200e3, 300e3]
    d = DensityProfile(h, [1e18, 1e17, 1e16]; origin = "narrow")
    e = ElectronProfile(h, [1e10, 1e11, 1e12], [200.0, 300.0, 400.0]; origin = "narrow")

    # Sampling within the native range is silent
    @test_nowarn d([150e3, 250e3])
    @test_nowarn e([150e3, 250e3])

    # Sampling past either end is extrapolation and must say so
    @test_logs (:warn, r"outside its native altitude range") d([50e3, 150e3])
    @test_logs (:warn, r"outside its native altitude range") e([150e3, 400e3])
end

@testitem "The model's neutrals argument must describe a whole atmosphere" begin
    electrons = ElectronProfile([90e3, 200e3, 400e3], [1e10, 1e11, 1e12],
                                [200.0, 500.0, 1500.0])

    # A single density source is rejected: it would otherwise silently become the density of
    # N2, O2, and O alike.
    @test_throws "neutrals must be a NeutralAtmosphere" AuroraModel(
        [100, 200], 180:-90:0, 100, @law(h -> fill(1e18, length(h))), electrons)

    # The swapped-arguments mistake keeps its dedicated hint
    @test_throws "Did you swap" AuroraModel([100, 200], 180:-90:0, 100, electrons, electrons)
end

@testitem "An electron source must be reproducible" begin
    z = make_altitude_grid(100, 200)
    iri_file = find_iri_file(; verbose = false)

    # A bare anonymous function cannot be saved to physics_state.jld2, so it is rejected at
    # construction, exactly as a species' density_source is.
    @test_throws ArgumentError Ionosphere(h -> (fill(1e11, length(h)),
                                                fill(200.0, length(h))), z)

    # A path and an ElectronProfile are both fine
    @test Ionosphere(iri_file, z) isa Ionosphere
    @test Ionosphere(read_iri_file(iri_file), z) isa Ionosphere
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

        # An export in other units must be rejected rather than silently rescaled, and the
        # error must name the unit mismatch — not imply the file is not a CCMC export.
        iri_m3 = joinpath(dir, "iri_m3.txt")
        write(iri_m3, "km Ne/m-3 Te/K\n100.0 4.8077E+10 159\n")
        @test_throws "no \"Ne/cm-3\" column" read_ccmc_iri(iri_m3)

        # A file with no Ne column at all is what "not a CCMC ModelWeb export" is reserved for
        @test_throws "Is this a CCMC ModelWeb export?" read_ccmc_iri(no_species)
    end
end

@testitem "read_ccmc_msis survives a short first data row" begin
    # The data column count must come from the modal token count across data rows, not just
    # the first one: a single short row after the header (e.g. a truncated first altitude)
    # must not make the reader think every row is that short and reject the whole file.
    mktempdir() do dir
        file = joinpath(dir, "short_first_row.txt")
        write(file, """
            Preamble line, not data
            Heit(km) N2den(cm-3) Oden(cm-3) O2den(cm-3)
            100.0    6.884E+12   4.793E+11
            150.0    2.690E+10   7.917E+09  1.965E+09
            200.0    3.127E+09   2.100E+09  1.338E+08
            """)
        p = read_ccmc_msis(file)
        @test p isa NeutralAtmosphere
        @test sort(collect(keys(p))) == [:N2, :O, :O2]
        # The short (100 km) row has too few tokens to parse any column and is dropped
        @test minimum(p[:N2].h) == 150e3
    end
end

@testitem "read_ccmc_msis survives CCMC's run-together MSIS00 header" begin
    # CCMC's NRLMSISE-00 export writes "Heden(cm-3)Arden(cm-3)" with no separating space, so
    # the header has one fewer token than the data has columns. Names before the glued token
    # still line up; names after it would return the neighbouring species' density.
    file = joinpath(@__DIR__, "test_data", "msis00_ccmc_output.txt")
    p = @test_warn "run-together name" read_ccmc_msis(file)

    # The species AURORA uses sit before the defect and are read correctly
    @test sort(collect(keys(p))) == [:N2, :O, :O2]
    @test p[:N2].n[1] ≈ 2.103e19 * 1e6 rtol=1e-12
    @test p[:O2].n[1] ≈ 5.641e18 * 1e6 rtol=1e-12
    @test p[:O].n[1]  ≈ 5.698e11 * 1e6 rtol=1e-12

    # Everything at or after the glued name is dropped rather than mislabelled. Without the
    # guard, :H would hold argon (9.697E+10) and :N would hold hydrogen (2.494E+07).
    for species in (:He, :Ar, :H, :N)
        @test !haskey(p, species)
    end

    # This export marks missing values with 0.000E+00 rather than 9.999E-38; O is absent
    # below 100 km either way
    @test minimum(p[:O].h) == 100e3
    @test minimum(p[:N2].h) == 0.0
end
