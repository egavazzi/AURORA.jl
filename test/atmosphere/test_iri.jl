# Tests for IRI file I/O, focusing on sentinel-value handling in load_iri_data.
#
# The reference file iri_20221102-1700_70N-19E.txt has sentinel -1 values in
# ne and Te at altitudes 60–79 km (20 bottom levels); all higher altitudes are valid.
#
# The reference file iri_20051008-2200_70N-19E.txt has no sentinel values and
# should load cleanly without warnings.

@testitem "load_iri_data: sentinel at bottom levels triggers warning and trims" begin
    iri_file = joinpath(@__DIR__, "test_data/", "iri_20221102-1700_70N-19E.txt")

    # A warning must be emitted
    @test_logs (:warn, r"sentinel -1") AURORA.load_iri_data(iri_file)

    data = AURORA.load_iri_data(iri_file) # (will throw a warning in REPL, nothing to worry)

    # The 20 bottom sentinel levels (60–79 km) must have been removed:
    # original grid is 60:1:700 (641 pts), valid start is 80 km → 621 pts remain.
    @test data.height_km[1] == 80.0
    @test length(data.height_km) == 621

    # No sentinel -1 values remain in ne or Te
    @test !any(data.ne .== -1)
    @test !any(data.Te .== -1)
end

@testitem "load_iri_data: clean file loads without warning" begin
    iri_file = joinpath(@__DIR__, "test_data", "iri_20051008-2200_70N-19E.txt")

    # No warning expected for a file without boundary sentinels
    @test_nowarn AURORA.load_iri_data(iri_file)
end

@testitem "run_iri: sentinel levels are trimmed, not stored" begin
    z = make_altitude_grid(80, 600)

    # Reaching into the D-region makes IRI return -1 for the lowest levels. They must be
    # dropped before they reach the log-space interpolation, exactly as when reading a file.
    p = @test_logs (:warn, r"sentinel -1") match_mode = :any begin
        run_iri(; year = 2005, month = 10, day = 8, hour = 22, minute = 0,
                  lat = 69.58, lon = 19.23, height = 50:5:700, verbose = false)
    end

    @test !any(p.ne .== -1)
    @test !any(p.Te .== -1)
    @test all(p.ne .> 0) && all(p.Te .> 0)
    @test p.h[1] > 50e3           # the sentinel levels at the bottom were removed

    ne, Te = p(z)
    @test all(isfinite, ne) && all(ne .> 0)
    @test all(isfinite, Te) && all(Te .> 0)
end
