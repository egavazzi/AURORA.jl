using DelimitedFiles: readdlm, writedlm
using Printf: @sprintf
using StyledStrings: @styled_str

############################################################################################
# Looking for files
############################################################################################
"""
    search_existing_iri_file(; year, month, day, hour, minute, lat, lon, height)

Search for an existing IRI data file matching the specified parameters.

This function scans the `internal_data/data_electron/` directory for IRI files
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
function search_existing_iri_file(; year, month, day, hour, minute, lat, lon, height,
                                   verbose=true)
    data_electron_directory = pkgdir(AURORA, "internal_data", "data_electron")
    data_electron_files = readdir(data_electron_directory)
    iri_files = data_electron_files[contains.(data_electron_files, "iri")]
    iri_files_fullpath = joinpath.(data_electron_directory, iri_files)
    verbose && print("Looking for an iri file that matches the parameters...")

    for (i, file) in enumerate(iri_files_fullpath)
        if iri_files[i][4] == '_' # Check that the file is not an old iri file
            year_file = tryparse(Int, iri_files[i][5:8])
            month_file = tryparse(Int, iri_files[i][9:10])
            day_file = tryparse(Int, iri_files[i][11:12])
            hour_file = tryparse(Int, iri_files[i][14:15])
            minute_file = tryparse(Int, iri_files[i][16:17])
            # First we do a pre-check, to avoid loading files unnecessarily
            if all([year_file, month_file, day_file, hour_file, minute_file] .==
                   [year, month, day, hour, minute])
                parameters = (; year, month, day, hour, minute, lat, lon, height)
                parameters_file = load_parameters_iri(file)
                # Now we check if all the parameters are the same
                if parameters == parameters_file
                    # If it is the case, that will be the file to load
                    verbose && println(" found one!")
                    verbose && println("   ∟ at $file")
                    return file
                end
            end
        end
    end

    verbose && println(" no file was found.")
    return nothing
end

############################################################################################
# Loading from file
############################################################################################
"""
    load_iri(iri_file)

Load an IRI file and return both parameters and data.

# Arguments
- `iri_file`: Path to the IRI file to load

# Returns
A NamedTuple with two fields:
- `parameters`: NamedTuple with (year, month, day, hour, minute, lat, lon, height)
- `data`: NamedTuple with all the IRI data columns
"""
function load_iri(iri_file)
    parameters = load_parameters_iri(iri_file)
    data = load_iri_data(iri_file)

    return (; parameters, data)
end

"""
    load_parameters_iri(iri_file)

Load calculation parameters from an IRI data file header.

Reads the header section of an IRI file and extracts the input parameters
that were used for the IRI model calculation.

# Arguments
- `iri_file::String`: Path to the IRI data file

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
function load_parameters_iri(iri_file)
    # Helper function to parse a value that could be Int or Float
    parse_numeric(s) = something(tryparse(Int64, s), parse(Float64, s))

    lines = open(iri_file, "r") do io
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
    load_iri_data(iri_file)

Load IRI model data from a file into a structured NamedTuple.

Reads an IRI data file and organizes all columns into a NamedTuple with
descriptive field names for easy access to densities, temperatures, and
other ionospheric parameters.

# Arguments
- `iri_file::String`: Path to the IRI data file

# Returns
- `NamedTuple`: IRI data with the following fields (all vectors except where noted):
  - `height_km`: Altitude grid (km)
  - `ne`: Electron density (m⁻³)
  - `Tn`: Neutral temperature (K)
  - `Ti`: Ion temperature (K)
  - `Te`: Electron temperature (K)
  - `nO⁺`, `nH⁺`, `nHe⁺`, `nO2⁺`, `nNO⁺`, `nCI`, `nN⁺`: Ion densities (m⁻³)
  - `NmF2`, `NmF1`, `NmE`: Peak densities at F2, F1, E layers
  - `hmF2`, `hmF1`, `hmE`: Peak heights at F2, F1, E layers (km)
  - `TEC`: Total Electron Content
  - `EqVertIonDrift`: Equatorial vertical ion drift
  - `foF2`: F2 critical frequency

# Notes
- All profile data are vectors with length equal to the iri altitude grid
"""
function load_iri_data(iri_file)
    # Find the header line with column names
    lines = readlines(iri_file)
    data_start_idx = findfirst(line -> startswith(line, "height(km)"), lines)

    # Use that information to load only the data
    data_matrix = readdlm(iri_file, skipstart = data_start_idx)

    # Extract columns and create named tuple
    # Column names from the header
    iri_data = (height_km = data_matrix[:, 1],
                ne = data_matrix[:, 2],          # electron density (m⁻³)
                Tn = data_matrix[:, 3],          # neutral temperature (K)
                Ti = data_matrix[:, 4],          # ion temperature (K)
                Te = data_matrix[:, 5],          # electron temperature (K)
                nO⁺ = data_matrix[:, 6],         # O⁺ density (m⁻³)
                nH⁺ = data_matrix[:, 7],         # H⁺ density (m⁻³)
                nHe⁺ = data_matrix[:, 8],        # He⁺ density (m⁻³)
                nO2⁺ = data_matrix[:, 9],        # O₂⁺ density (m⁻³)
                nNO⁺ = data_matrix[:, 10],       # NO⁺ density (m⁻³)
                nCI = data_matrix[:, 11],        # Cluster ions density (m⁻³)
                nN⁺ = data_matrix[:, 12],        # N⁺ density (m⁻³)
                NmF2 = data_matrix[:, 13],       # F2 peak density
                hmF2 = data_matrix[:, 14],       # F2 peak height
                NmF1 = data_matrix[:, 15],       # F1 peak density
                hmF1 = data_matrix[:, 16],       # F1 peak height
                NmE = data_matrix[:, 17],        # E peak density
                hmE = data_matrix[:, 18],        # E peak height
                TEC = data_matrix[:, 19],        # Total Electron Content
                EqVertIonDrift = data_matrix[:, 20],  # Equatorial vertical ion drift
                foF2 = data_matrix[:, 21])

    # Validate: check that the loaded data doesn't contain only -1 sentinel values
    if all(iri_data.ne .== -1) && all(iri_data.Te .== -1)
        error("The IRI file at\n  $(iri_file)\n" *
              "contains only -1 sentinel values (no valid ionospheric profiles).\n" *
              "You might want to use another file or regenerate it.")
    end

    # Check for sentinel -1 values at the lowest and highest altitudes.
    # If found, emit a warning, drop those and let the interpolator fill them.
    sentinel_mask = (iri_data.ne .== -1) .| (iri_data.Te .== -1)
    if any(sentinel_mask)
        first_valid = findfirst(!, sentinel_mask)
        last_valid  = findlast(!, sentinel_mask)
        n_bottom = first_valid - 1
        n_top    = length(sentinel_mask) - last_valid

        if n_bottom > 0 || n_top > 0
            h = iri_data.height_km
            if n_bottom > 0 && n_top > 0
                location_str = "the $(n_bottom) lowest altitude levels " *
                               "($(h[1]) – $(h[n_bottom]) km) and " *
                               "$(n_top) highest altitude levels " *
                               "($(h[last_valid + 1]) – $(h[end]) km)"
            elseif n_bottom > 0
                location_str = "the $(n_bottom) lowest altitude level(s) " *
                               "($(h[1]) – $(h[n_bottom]) km)"
            else
                location_str = "the $(n_top) highest altitude level(s) " *
                               "($(h[last_valid + 1]) – $(h[end]) km)"
            end
            @warn "IRI file at\n  $(iri_file)\n" *
                  "has sentinel -1 values in ne or Te at $(location_str).\n" *
                  "These levels will be dropped and the interpolator will fill them " *
                  "through extrapolation."

            valid_range = first_valid:last_valid
            iri_data = map(v -> v isa AbstractVector ? v[valid_range] : v, iri_data)
        end
    end

    return iri_data
end

############################################################################################
# Saving to file
############################################################################################
"""
    save_iri_data(iri_data, parameters)

Save IRI model data to a text file with metadata header.

Creates a formatted text file containing the IRI calculation parameters in the header
followed by the data matrix. The filename is automatically generated based on the
input parameters, and if a file with the same name exists, a unique name is created.

# Arguments
- `iri_data::Matrix`: IRI data matrix from `calculate_iri_data()` (with header row)
- `parameters::NamedTuple`: Parameters used for IRI calculation, must contain:
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
`iri_YYYYMMDD-HHMM_LATN-LONE.txt`

# Notes
- Files are saved to `internal_data/data_electron/` directory
- Existing files are not overwritten; a suffix is added to the filename
"""
function save_iri_data(iri_data, parameters; verbose=true)
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
    filename = "iri_$year$month_str$day_str-$(hour_str)$(minute_str)_$(lat)N-$(lon)E.txt"
    directory = pkgdir(AURORA, "internal_data", "data_electron")
    fullpath = joinpath(directory, filename)
    fullpath = rename_if_exists(fullpath) # to avoid writing over files
    filename = splitpath(fullpath)[end]
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
        writedlm(f, iri_data)
    end
    verbose && println(styled"File {bold:$filename} created under {underline:$directory}.")

    return fullpath
end
