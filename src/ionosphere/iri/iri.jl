include("calculation.jl")
include("io.jl")
include("electron_profile.jl")

"""
    find_iri_file(; year=2018, month=12, day=7, hour=11, minute=15,
                    lat=76, lon=5, height=85:1:700)

Find or create an IRI model data file for the specified conditions.

This is the cached, file-based route to the IRI-2020 model. It first searches
`internal_data/data_electron/` for an existing IRI file matching the given parameters. If no
matching file is found, it calculates new IRI data using the Python iri2020 package and saves
it there, so that a later call with the same parameters reuses the file instead of running the
model again. The iri2020 package will compile and run some fortran code under the hood.

`read_iri_file(find_iri_file(; ...))` turns the result into an [`ElectronProfile`](@ref).
[`run_iri`](@ref) runs the model directly and returns that profile, writing a file only if
asked to.

# Keyword Arguments
- `year::Int=2018`: Year
- `month::Int=12`: Month (1-12)
- `day::Int=7`: Day of month (1-31)
- `hour::Int=11`: Hour in Universal Time (0-23)
- `minute::Int=15`: Minute (0-59)
- `lat::Real=76`: Geographic latitude in degrees North
- `lon::Real=5`: Geographic longitude in degrees East
- `height::AbstractRange=85:1:700`: Altitude range in km

# Returns
- `String`: Full path to the IRI data file

# Notes
- Default parameters correspond to the VISIONS-2 rocket launch conditions
- Files are stored in `internal_data/data_electron/` directory
"""
function find_iri_file(;
                       year = 2018,
                       month = 12,
                       day = 7,
                       hour = 11,
                       minute = 15,
                       lat = 76,
                       lon = 5,
                       height = 85:1:700,
                       verbose = true)

    # First check if we have an iri file with these parameters
    file_to_load = search_existing_iri_file(; year, month, day, hour, minute, lat, lon, height,
                                            verbose)
    if !isnothing(file_to_load)
        return file_to_load
    end

    # Otherwise, calculate new IRI data and save it, so that the next call with these
    # parameters finds it.
    iri_data, parameters = calculate_iri_data(; year, month, day, hour, minute, lat, lon,
                                              height, verbose)
    file_to_load = save_iri_data(iri_data, parameters; verbose)

    return file_to_load
end
