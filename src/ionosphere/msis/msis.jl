include("calculation.jl")
include("io.jl")

"""
    find_msis_file(; year=2018, month=12, day=7, hour=11, minute=15,
                    lat=76, lon=5, height=85:1:700)

Find or create a MSIS model data file for the specified conditions.

This is the cached, file-based route to the NRLMSIS 2.1 model. It first searches
`internal_data/data_neutrals/` for an existing MSIS file matching the given parameters. If no
matching file is found, it calculates new MSIS data using the Python pymsis package and saves
it there, so that a later call with the same parameters reuses the file instead of running the
model again. The pymsis package will download, compile and run some fortran code under the
hood.

`read_msis_file(find_msis_file(; ...))` turns the result into a [`NeutralAtmosphere`](@ref).
[`run_msis`](@ref) runs the model directly and returns that atmosphere, writing a file only if
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
- `String`: Full path to the MSIS data file

# Notes
- Default parameters correspond to the VISIONS-2 rocket launch conditions
- Files are stored in `internal_data/data_neutrals/` directory
"""
function find_msis_file(;
                        year = 2018,
                        month = 12,
                        day = 7,
                        hour = 11,
                        minute = 15,
                        lat = 76,
                        lon = 5,
                        height = 85:1:700,
                        verbose = true)

    # First check if we have a msis file with these parameters
    file_to_load = search_existing_msis_file(; year, month, day, hour, minute, lat, lon, height,
                                             verbose)
    if !isnothing(file_to_load)
        return file_to_load
    end

    # Otherwise, calculate new MSIS data and save it, so that the next call with these
    # parameters finds it.
    msis_data, parameters = calculate_msis_data(; year, month, day, hour, minute, lat, lon,
                                                height, verbose)
    file_to_load = save_msis_data(msis_data, parameters; verbose)

    return file_to_load
end

# Deprecated alias for backward compatibility
function find_nrlmsis_file(;
                           year = 2018,
                           month = 12,
                           day = 7,
                           hour = 11,
                           minute = 15,
                           lat = 76,
                           lon = 5,
                           height = 85:1:700,
                           verbose = true)
    @warn "find_nrlmsis_file() is deprecated, use find_msis_file() instead" maxlog = 1
    return find_msis_file(; year, month, day, hour, minute, lat, lon, height, verbose)
end
