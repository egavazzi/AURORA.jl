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
    @test_warn r"sentinel -1" AURORA.load_iri_data(iri_file)

    data = AURORA.load_iri_data(iri_file)

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
