using DelimitedFiles: readdlm, writedlm
using Printf: @sprintf
using Term: @bold, @underline

############################################################################################
# Looking for files
############################################################################################
"""
    search_existing_msis_file(; year, month, day, hour, minute, lat, lon, height)

Search for an existing MSIS data file matching the specified parameters.

This function scans the `internal_data/data_neutrals/` directory for MSIS files
with matching time, location, and altitude grid parameters. It performs a quick
pre-check on filenames before loading and comparing full parameters.

# Keyword Arguments
- `year::Int`: Year
- `month::Int`: Month (1-12)
- `day::Int`: Day of month (1-31)
- `hour::Int`: Hour in Universal Time (0-23)
- `minute::Int`: Minute (0-59)
- `lat::Real`: Geographic latitude in degrees North
- `lon::Real`: Geographic longitude in degrees East
- `height::AbstractRange`: Altitude range in km

# Returns
- `Union{String, Nothing}`: Full path to matching file, or `nothing` if not found
"""
function search_existing_msis_file(; year, month, day, hour, minute, lat, lon, height)
    data_neutrals_directory = pkgdir(AURORA, "internal_data", "data_neutrals")
    data_neutrals_files = readdir(data_neutrals_directory)
    msis_files = data_neutrals_files[contains.(data_neutrals_files, "msis")]
    msis_files_fullpath = joinpath.(data_neutrals_directory, msis_files)
    print("Looking for an msis file that matches the parameters...")

    for (i, file) in enumerate(msis_files_fullpath)
        if msis_files[i][5] == '_' # Check that the file is not an old msis file
            year_file = tryparse(Int, msis_files[i][6:9])
            month_file = tryparse(Int, msis_files[i][10:11])
            day_file = tryparse(Int, msis_files[i][12:13])
            hour_file = tryparse(Int, msis_files[i][15:16])
            minute_file = tryparse(Int, msis_files[i][17:18])
            # First we do a pre-check, to avoid loading files unnecessarily
            if all([year_file, month_file, day_file, hour_file, minute_file] .==
                   [year, month, day, hour, minute])
                parameters = (; year, month, day, hour, minute, lat, lon, height)
                parameters_file = load_parameters_msis(file)
                # Now we check if all the parameters are the same
                if parameters == parameters_file
                    # If it is the case, that will be the file to load
                    println(" found one!")
                    println("   ∟ at $file")
                    return file
                end
            end
        end
    end

    println(" no file was found.")
    return nothing
end

############################################################################################
# Loading from file
############################################################################################
"""
    load_msis(msis_file)

Load a MSIS file and return both parameters and data.

# Arguments
- `msis_file`: Path to the MSIS file to load

# Returns
A NamedTuple with two fields:
- `parameters`: NamedTuple with (year, month, day, hour, minute, lat, lon, height)
- `data`: NamedTuple with all the MSIS data columns
"""
function load_msis(msis_file)
    parameters = load_parameters_msis(msis_file)
    data = load_msis_data(msis_file)

    return (; parameters, data)
end

"""
    load_parameters_msis(msis_file)

Load calculation parameters from a MSIS data file header.

Reads the header section of a MSIS file and extracts the input parameters
that were used for the MSIS model calculation.

# Arguments
- `msis_file::String`: Path to the MSIS data file

# Returns
- `NamedTuple`: Parameters with fields:
  - `year::Int`: Year
  - `month::Int`: Month (1-12)
  - `day::Int`: Day (1-31)
  - `hour::Int`: Hour (0-23)
  - `minute::Int`: Minute (0-59)
  - `lat::Real`: Latitude (degrees North)
  - `lon::Real`: Longitude (degrees East)
  - `height::AbstractRange`: Altitude range (km)
"""
function load_parameters_msis(msis_file)
    # Helper function to parse a value that could be Int or Float
    parse_numeric(s) = something(tryparse(Int64, s), parse(Float64, s))

    lines = open(msis_file, "r") do io
        [readline(io) for _ in 1:12]
    end

    # Parse parameters from the header
    # Expected format: "parameter_name = value"
    year = parse(Int, split(lines[3], "=")[2])
    month = parse(Int, split(lines[4], "=")[2])
    day = parse(Int, split(lines[5], "=")[2])
    hour = parse(Int, split(lines[6], "=")[2])
    minute = parse(Int, split(lines[7], "=")[2])
    lat = parse_numeric(strip(split(lines[8], "=")[2]))
    lon = parse_numeric(strip(split(lines[9], "=")[2]))

    # Parse the height range (format: "start:step:stop")
    height_str = strip(split(lines[10], "=")[2])
    height_parts = split(height_str, ":")
    height_start = parse_numeric(strip(height_parts[1]))
    height_step = parse_numeric(strip(height_parts[2]))
    height_stop = parse_numeric(strip(height_parts[3]))
    height = height_start:height_step:height_stop

    return (; year, month, day, hour, minute, lat, lon, height)
end

"""
    load_msis_data(msis_file)

Load MSIS model data from a file into a structured NamedTuple.

Reads a MSIS data file and organizes all columns into a NamedTuple with
descriptive field names for easy access to densities, temperature, and
other atmospheric parameters.

# Arguments
- `msis_file::String`: Path to the MSIS data file

# Returns
- `NamedTuple`: MSIS data with the following fields (all vectors):
  - `height_km`: Altitude grid (km)
  - `air`: Total mass density (kg/m³)
  - `N2`: Molecular nitrogen density (m⁻³)
  - `O2`: Molecular oxygen density (m⁻³)
  - `O`: Atomic oxygen density (m⁻³)
  - `He`: Helium density (m⁻³)
  - `H`: Atomic hydrogen density (m⁻³)
  - `Ar`: Argon density (m⁻³)
  - `N`: Atomic nitrogen density (m⁻³)
  - `anomalousO`: Anomalous oxygen density (m⁻³)
  - `NO`: Nitric oxide density (m⁻³)
  - `T`: Temperature (K)

# Notes
- All data are vectors with length equal to the msis altitude grid
"""
function load_msis_data(msis_file)
    # Find the header line with column names
    lines = readlines(msis_file)
    data_start_idx = findfirst(line -> startswith(line, "height(km)"), lines)

    # Use that information to load only the data
    data_matrix = readdlm(msis_file, skipstart = data_start_idx)

    # Extract columns and create named tuple
    # Column names from the header
    msis_data = (height_km = data_matrix[:, 1],
                 air = data_matrix[:, 2],           # total mass density (kg/m³)
                 N2 = data_matrix[:, 3],            # N₂ density (m⁻³)
                 O2 = data_matrix[:, 4],            # O₂ density (m⁻³)
                 O = data_matrix[:, 5],             # O density (m⁻³)
                 He = data_matrix[:, 6],            # He density (m⁻³)
                 H = data_matrix[:, 7],             # H density (m⁻³)
                 Ar = data_matrix[:, 8],            # Ar density (m⁻³)
                 N = data_matrix[:, 9],             # N density (m⁻³)
                 anomalousO = data_matrix[:, 10],   # anomalous O density (m⁻³)
                 NO = data_matrix[:, 11],           # NO density (m⁻³)
                 T = data_matrix[:, 12])            # Temperature (K)

    return msis_data
end

############################################################################################
# Saving to file
############################################################################################
"""
    save_msis_data(msis_data, parameters)

Save MSIS model data to a text file with metadata header.

Creates a formatted text file containing the MSIS calculation parameters in the header
followed by the data matrix. The filename is automatically generated based on the
input parameters, and if a file with the same name exists, a unique name is created.

# Arguments
- `msis_data::Matrix`: MSIS data matrix from `calculate_msis_data()` (with header row)
- `parameters::NamedTuple`: Parameters used for MSIS calculation, must contain:
  - `year`, `month`, `day`, `hour`, `minute`: Time specification
  - `lat`, `lon`: Location (degrees)
  - `height`: Altitude range (km)

# Returns
- `String`: Full path to the created file

# File Format
The file contains:
1. Header section with input parameters
2. Column headers
3. Data matrix (one row per altitude)

# Filename Convention
`msis_YYYYMMDD-HHMM_LATN-LONE.txt`

# Notes
- Files are saved to `internal_data/data_neutrals/` directory
- Existing files are not overwritten; a suffix is added to the filename
"""
function save_msis_data(msis_data, parameters)
    # Unpack the parameters
    year = parameters.year
    month = parameters.month
    day = parameters.day
    hour = parameters.hour
    minute = parameters.minute
    lat = parameters.lat
    lon = parameters.lon
    height = parameters.height
    # Make strings with 2 digits (for the file name)
    month_str = @sprintf "%02d" month
    day_str = @sprintf "%02d" day
    hour_str = @sprintf "%02d" hour
    minute_str = @sprintf "%02d" minute

    # Make filename
    filename = "msis_$year$month_str$day_str-$(hour_str)$(minute_str)_$(lat)N-$(lon)E.txt"
    directory = pkgdir(AURORA, "internal_data", "data_neutrals")
    fullpath = joinpath(directory, filename)
    fullpath = rename_if_exists(fullpath) # to avoid writing over files
    filename = splitpath(fullpath)[end] # update filename as fullpath has been updated
    # Write to the file
    open(fullpath, "w") do f
        write(f, "Input parameters:\n")
        write(f, "\n")
        write(f, "year = $year\n")
        write(f, "month = $month\n")
        write(f, "day = $day\n")
        write(f, "hour = $hour\n")
        write(f, "minute = $minute\n")
        write(f, "latitude = $lat\n")
        write(f, "longitude = $lon\n")
        write(f, "height (km) = $height\n")
        write(f, "time_type = Universal Time\n")
        write(f, "coordinate_type = Geographic\n")
        write(f, "\n")
        writedlm(f, msis_data)
    end
    println("File " * @bold("$filename") * " created under " * @underline("$directory") *
            ".")

    return fullpath
end
