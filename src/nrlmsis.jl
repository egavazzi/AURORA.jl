using HTTP
using JSON
using Downloads
using DelimitedFiles
using ProgressMeter
using Dates
using Term
using Printf

## TODO : add docstrings

# Defaults to Visions2 launch conditions
function find_nrlmsis_file(;
    year = 2018,
    month = 12,
    day = 7,
    hour = 11,
    minute = 15,
    lat = 76,
    lon = 5,
    height = 85:1:700,
    )

    # First check if we have an msis file with these parameters
    data_neutrals_directory = pkgdir(AURORA, "internal_data", "data_neutrals")
    data_neutrals_files = readdir(data_neutrals_directory)
    msis_files = data_neutrals_files[contains.(data_neutrals_files, "msis")]
    msis_files_fullpath = joinpath.(data_neutrals_directory, msis_files)
    print("Looking for an msis file that matches the parameters...")
    found_it = 0
    file_to_load = ""
    for (i, file) in enumerate(msis_files_fullpath)
        if msis_files[i][5] == '_' # Check that the file is not an old msis file
            year_file = tryparse(Int, msis_files[i][6:9])
            month_file = tryparse(Int, msis_files[i][10:11])
            day_file = tryparse(Int, msis_files[i][12:13])
            hour_file = tryparse(Int, msis_files[i][15:16])
            minute_file = tryparse(Int, msis_files[i][17:18])
            # First we do a pre-check, to avoid loading files unnecessarily
            if all([year_file, month_file, day_file, hour_file, minute_file] .== [year, month, day, hour, minute])
                parameters = (; year, month, day, hour, minute, lat, lon, height)
                parameters_file = load_parameters_msis(file)
                # Now we check if all the parameters are the same
                if parameters == parameters_file
                    # If it is the case, that will be the file to load
                    found_it = 1
                    file_to_load = file
                    println(" found one!")
                end
            end
        end
    end

    # If we did not find a file matching the parameters, we need to calculate one
    if found_it == 0
        println(" no file was found.")
        # nrlmsis_data, parameters = download_msis_data(year, month, day, hour, minute, lat, lon, height)
        nrlmsis_data, parameters = calculate_msis_data(year, month, day, hour, minute, lat, lon, height)
        # and save it
        file_to_load = save_msis_data(nrlmsis_data, parameters)
        found_it = 1
    end

    if found_it == 1
        return file_to_load
    else
        error("This is not normal, something is wrong with find_nrlmsis_file()")
    end
end



function load_parameters_msis(nrlmsis_file)
    local data
    # Loading the `parameters` section of the file
    open(`head -n12 $nrlmsis_file`) do io
        data = readdlm(io)
    end
    # Extracting the parameters
    year = data[2, 3]
    month = data[3, 3]
    day = data[4, 3]
    hour = data[5, 3]
    minute = data[6, 3]
    lat = data[7, 3]
    lon = data[8, 3]
    height_str = data[9, 4]
    height_start = split(height_str, ":")[1]
    height_step = split(height_str, ":")[2]
    height_stop = split(height_str, ":")[3]
    # height is loaded as a string
    # We try to convert to Int64 and if it fails we try to Float64
    try
        height_start = parse(Int64, height_start)
        height_step = parse(Int64, height_step)
        height_stop = parse(Int64, height_stop)
    catch
        height_start = parse(Float64, height_start)
        height_step = parse(Float64, height_step)
        height_stop = parse(Float64, height_stop)
    end
    height = height_start:height_step:height_stop

    parameters_file = (; year, month, day, hour, minute, lat, lon, height)
    return parameters_file
end


# This method is depreciated. The HTTP request frequently has issues with the server.
function download_msis_data(year = 2018, month = 12, day = 7, hour = 11, minute = 15,
    lat = 76, lon = 5, height = 85:1:700)

    println("Starting to download msis data from https://kauai.ccmc.gsfc.nasa.gov/instantrun/nrlmsis/. It can take a few minutes (on the server side)...")

    # Initialize matrix
    nrlmsis_data = Matrix{any}

    # URL to the nrlmsis api
    url = "https://kauai.ccmc.gsfc.nasa.gov/instantrun/api/nrlmsis"

    # HEADERS, should not be changed
    headers = Dict(
        "Host" => "kauai.ccmc.gsfc.nasa.gov",
        "User-Agent" => "Mozilla/5.0 (X11; Linux x86_64; rv:109.0) Gecko/20100101 Firefox/119.0",
        "Accept" => "application/json",
        "Accept-Language" => "no,sv;q=0.8,fr;q=0.6,en-GB;q=0.4,en;q=0.2",
        "Accept-Encoding" => "gzip, deflate, br",
        "Referer" => "https://kauai.ccmc.gsfc.nasa.gov/instantrun/nrlmsis/",
        # "Content-Length" => "229",
        "Origin" => "https://kauai.ccmc.gsfc.nasa.gov",
        "DNT" => "1",
        "Connection" => "keep-alive",
        "Sec-Fetch-Dest" => "empty",
        "Sec-Fetch-Mode" => "no-cors",
        "Sec-Fetch-Site" => "same-origin",
        "content-type" => "application/json",
        "Pragma" => "no-cache",
        "Cache-Control" => "no-cache",
    )
    # BODY, that's the part that changes
    # format the time
    datetime = DateTime(year, month, day, hour, minute)
    datetime = Dates.format(datetime, "yyyy-mm-ddTHH:MM:SS.000Z")

    # The new nrlmsis api does not support height ranges with a length greater than 149. So
    # if length(height) > 149, we need to do several quary and concatenate the results.
    n_query = ceil(Int, length(height) / 149) # number of query necessary
    for i_query in 1:n_query
        if i_query < n_query
            local_height = height[(1:149) .+ (i_query - 1) * 149]
        elseif i_query == n_query
            local_height = height[(i_query - 1) * 149 + 1:end]
        end
        body = Dict(
            "version" => "NRLMSIS2",
            "timeType" => "ut",         # universal time
            "coordinateType" => "geo",  # geographic
            "datetime" => datetime,
            "lat" => lat,
            "lon" => lon,
            "height" => 100,
            "profileType" => 1,       # height profile
            "start" => local_height[1],
            "stop" => local_height[end],
            "step" => step(local_height),
            "fdaily" => -1,           # default (observational data)
            "f3month" => -1,          # default (observational data)
            "ap" => -1,               # default (observational data)
            "_datelimit" => true,
        )
        body = JSON.json(body) # convert to JSON
        output_file_link = HTTP.request("POST", url; headers, body, connect_timeout = 300)     # request results
        output_file_link = JSON.parse(String(output_file_link.body))    # from JSON to Julia
        output_file_link = output_file_link["txt"]                      # get only the link


        # read the results and put them in the `nrlmsis` matrix
        if i_query == 1
            # first one, we keep the headers (skipstart = 1)
            nrlmsis_data = readdlm(Downloads.download(output_file_link))
        else
            # we skip the headers (skipstart = 2)
            data = readdlm(Downloads.download(output_file_link), skipstart = 1)
            nrlmsis_data = vcat(nrlmsis_data, data)
        end
        println("Downloading the nrlmsis data ($i_query/$n_query).")
        # It seems like the nrlmsis server get capacity problems when we do multiple requests
        # following each other in a short time. So the idea here is to wait two seconds between
        # each request
        sleep(2)
    end

    parameters = (; year, month, day, hour, minute, lat, lon, height)
    return nrlmsis_data, parameters
end


using HTTP
using PythonCall
function calculate_msis_data(year = 2018, month = 12, day = 7, hour = 11, minute = 15,
    lat = 76, lon = 5, height = 85:1:700)

    print("Calculating msis data...")

    datetime = pyimport("datetime")
    time = datetime.datetime(year, month, day, hour, minute, 0)



    # import nrlmsis 2.1 model from the Python package 'pymsis'
    msis = pyimport("pymsis.msis")

    # Check if the SW-ALL.csv file is already downloaded or not. If it is not downloaded,
    # download it. This is normally done by the pymsis package, but due to some issues with
    # HTTP request and SSL certificates when calling the package from Julia, we do this
    # action directly in Julia instead.
    os = pyimport("os")
    path_to_pymsis = os.path.dirname(msis.__file__)
    path_to_pymsis = pyconvert(String, path_to_pymsis)
    if isfile(joinpath(path_to_pymsis, "SW-All.csv")) == false
        print(" downloading ap and F10.7 data from (https://celestrak.org/SpaceData/SW-All.csv)...")
        HTTP.download("https://celestrak.org/SpaceData/SW-All.csv", path_to_pymsis)
    end

    # run the model
    nrlmsis_data = msis.run(time, Py(lon), Py(lat), Py(height), geomagnetic_activity=Py(-1))
    # convert from Python array to Julia array
    nrlmsis_data = pyconvert(Array, nrlmsis_data) # array of size (1, 1, 1, n_z, 11)
    nrlmsis_data = dropdims(nrlmsis_data; dims = (1, 2, 3)) # convert to size (n_z, 11)
    # add a column with the altitude
    nrlmsis_data = hcat(Vector(height), nrlmsis_data)
    # add a header with the name of columns
    nrlmsis_data = vcat(["height(km)" "air(kg/m3)" "N2(m-3)" "O2(m-3)" "O(m-3)" "He(m-3)" "H(m-3)" "Ar(m-3)" "N(m-3)" "anomalousO(m-3)" "NO(m-3)" "T(K)"], nrlmsis_data)

    println(" done.")

    parameters = (; year, month, day, hour, minute, lat, lon, height)
    return nrlmsis_data, parameters
end


function save_msis_data(nrlmsis_data, parameters)
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
        # First we write the parameters
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
        # Then we write the data
        writedlm(f, nrlmsis_data)
    end
    println("File " * @bold("$filename") * " created under " * @underline("$directory") * ".")

    return fullpath
end





# m = nrlmsis_data[2, 2]
# d = nrlmsis_data[2, 3]
# hourdotminute = nrlmsis_data[2, 5]
# h = match(r"(.*)(?=\.)", string(hourdotminute)).match; h = parse(Int, h) # Could use split...
# m = match(r"(?<=\.)(.*)", string(hourdotminute)).match; m = parse(Int, m)
# lat = nrlmsis_data[2, 7]
# lon = nrlmsis_data[2, 8]
# height_start = nrlmsis_data[2, 6]
# height_end = nrlmsis_data[end, 6]
# height_step = nrlmsis_data[3, 6] - nrlmsis_data[2, 6]
# height = height_start:height_step:height_end
