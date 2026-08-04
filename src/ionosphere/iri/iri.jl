include("calculation.jl")
include("io.jl")
include("electron_profile.jl")

"""
    find_iri_file(; year=2018, month=12, day=7, hour=11, minute=15,
                    lat=76, lon=5, height=85:1:700)

Find or create an IRI model data file for the specified conditions.

It first searches for an existing IRI file matching the given parameters.
If no matching file is found, it calculates new IRI data using the Python iri2020 package
and saves it to a file. The iri2020 package will compile and run some fortran code under
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
                       height = 85:1:700,
                       verbose = true)

    # First check if we have an iri file with these parameters
    file_to_load = search_existing_iri_file(; year, month, day, hour, minute, lat, lon, height,
                                            verbose)
    if !isnothing(file_to_load)
        return file_to_load
    end

    # Otherwise, calculate and save new IRI data. Generating files is deprecated: prefer
    # run_iri, which returns an ElectronProfile stored in (and round-tripped with) the model
    # rather than written to disk.
    @warn """
    find_iri_file() is generating a new IRI file via the Python `iri2020` package, which is \
    deprecated. Prefer run_iri for new runs, e.g.:
        iri   = run_iri(; year=$year, month=$month, day=$day, hour=$hour, minute=$minute, lat=$lat, lon=$lon)
        model = AuroraModel(altitude_lims, θ_lims, E_max, neutrals, electrons)
    Reading existing IRI files (via read_iri_file / AuroraModel) remains supported.""" maxlog = 1
    iri_data, parameters = calculate_iri_data(; year, month, day, hour, minute, lat, lon,
                                              height, verbose)
    file_to_load = save_iri_data(iri_data, parameters; verbose)

    return file_to_load
end
