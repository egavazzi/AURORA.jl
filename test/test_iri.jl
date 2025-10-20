@testitem "Search for file that does not exist" begin
    PARAMETERS = [2009, 12, 10, 4, 5, 76, 5, 82:1:600]
    file_to_load = AURORA.search_existing_iri_file(;
                                            year = PARAMETERS[1],
                                            month = PARAMETERS[2],
                                            day = PARAMETERS[3],
                                            hour = PARAMETERS[4],
                                            minute = PARAMETERS[5],
                                            lat = PARAMETERS[6],
                                            lon = PARAMETERS[7],
                                            height = PARAMETERS[8])
    # This might fail on local if a file with these parameters has been previously saved
    @test isnothing(file_to_load)
end

@testitem "Search for file that already exists" begin
    PARAMETERS = [2009, 12, 10, 4, 5, 76, 5, 82:1:600]
    iri_data, iri_parameters = AURORA.calculate_iri_data(;
                                                  year = PARAMETERS[1],
                                                  month = PARAMETERS[2],
                                                  day = PARAMETERS[3],
                                                  hour = PARAMETERS[4],
                                                  minute = PARAMETERS[5],
                                                  lat = PARAMETERS[6],
                                                  lon = PARAMETERS[7],
                                                  height = PARAMETERS[8])
    AURORA.save_iri_data(iri_data, iri_parameters)
    file_to_load = AURORA.search_existing_iri_file(;
                                            year = PARAMETERS[1],
                                            month = PARAMETERS[2],
                                            day = PARAMETERS[3],
                                            hour = PARAMETERS[4],
                                            minute = PARAMETERS[5],
                                            lat = PARAMETERS[6],
                                            lon = PARAMETERS[7],
                                            height = PARAMETERS[8])
    @test !isnothing(file_to_load)
end

@testitem "Loading" begin
    iri_file = AURORA.find_iri_file()
    AURORA.load_iri(iri_file)
    @test true
end

@testitem "Interpolating" begin
    iri_file = find_iri_file()
    iri = AURORA.load_iri(iri_file)
    h_atm = make_altitude_grid(60, 500)
    # simple interpolation
    AURORA.interpolate_profile(iri.ne, iri.z_iri, h_atm; log_interpolation = true)
    AURORA.interpolate_profile(iri.Te, iri.z_iri, h_atm; log_interpolation = false)
    @test true
    # interpolate all iri data
    AURORA.interpolate_iri_to_grid(iri_data, h_atm)
    @test true
end
