include("calculation.jl")
include("io.jl")
include("interpolation.jl")

"""
    find_iri_file(; year=2018, month=12, day=7, hour=11, minute=15,
                    lat=76, lon=5, height=85:1:700)

Find or create an IRI model data file for the specified conditions.

It first searches for an existing IRI file matching the given parameters.
If no matching file is found, it calculates new IRI data using the Python iri2016 package
and saves it to a file. The iri2016 package will compile and run some fortran code under
the hood.

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
                       height = 85:1:700)

    # First check if we have an iri file with these parameters
    file_to_load = search_existing_iri_file(; year, month, day, hour, minute, lat, lon, height)
    if !isnothing(file_to_load)
        return file_to_load
    end

    # Otherwise, calculate and save new IRI data
    iri_data, parameters = calculate_iri_data(; year, month, day, hour, minute, lat, lon, height)
    file_to_load = save_iri_data(iri_data, parameters)

    return file_to_load
end
