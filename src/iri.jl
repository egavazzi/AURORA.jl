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
function find_iri_file(;
    year = 2018,
    month = 12,
    day = 7,
    hour = 11,
    minute = 15,
    lat = 76,
    lon = 5,
    height = 85:1:700,
    )

    # First check if we have an iri file with these parameters
    data_electron_directory = pkgdir(AURORA, "internal_data", "data_electron")
    data_electron_files = readdir(data_electron_directory)
    iri_files = data_electron_files[contains.(data_electron_files, "iri")]
    iri_files_fullpath = joinpath.(data_electron_directory, iri_files)
    print("Looking for an iri file that matches the parameters...")
    found_it = 0
    file_to_load = ""
    for (i, file) in enumerate(iri_files_fullpath)
        if iri_files[i][4] == '_' # Check that the file is not an old iri file
            year_file = tryparse(Int, iri_files[i][5:8])
            month_file = tryparse(Int, iri_files[i][9:10])
            day_file = tryparse(Int, iri_files[i][11:12])
            hour_file = tryparse(Int, iri_files[i][14:15])
            minute_file = tryparse(Int, iri_files[i][16:17])
            # First we do a pre-check, to avoid loading files unnecessarily
            if all([year_file, month_file, day_file, hour_file, minute_file] .== [year, month, day, hour, minute])
                parameters = (; year, month, day, hour, minute, lat, lon, height)
                parameters_file = load_parameters_iri(file)
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

    # If we did not find a file matching the parameters, we need to download one
    if found_it == 0
        println(" no file was found.")
        # iri_data, parameters = download_iri_data(year, month, day, hour, minute, lat, lon, height)
        iri_data, parameters = calculate_iri_data(year, month, day, hour, minute, lat, lon, height)
        # and save it
        file_to_load = save_iri_data(iri_data, parameters)
        found_it = 1
    end

    if found_it == 1
        return file_to_load
    else
        error("This is not normal, something is wrong with find_iri_file()")
    end
end



function download_iri_data(year=2018, month=12, day=7, hour=11, minute=15,
    lat = 76, lon = 5, height = 85:1:700)

    println("Starting to download iri data from https://kauai.ccmc.gsfc.nasa.gov/instantrun/iri/. It can take a few minutes (on the server side)...")

    # URL to the iri api
    url = "https://kauai.ccmc.gsfc.nasa.gov/instantrun/api/iri/2020/"

    # HEADERS, should not be changed
    headers = Dict(
        "Accept" => "application/json",
        "Accept-Encoding" => "gzip, deflate, br",
        "Accept-Language" => "no,sv;q=0.8,fr;q=0.6,en-GB;q=0.4,en;q=0.2",
        "Cache-Control" => "no-cache",
        "Connection" => "keep-alive",
        "Content-Type" => "application/json",
        "DNT" => "1",
        "Host" => "kauai.ccmc.gsfc.nasa.gov",
        "Origin" => "https://kauai.ccmc.gsfc.nasa.gov",
        "Pragma" => "no-cache",
        "Referer" => "https://kauai.ccmc.gsfc.nasa.gov/instantrun/iri/",
        "Sec-Fetch-Dest" => "empty",
        "Sec-Fetch-Mode" => "no-cors",
        "Sec-Fetch-Site" => "same-origin",
        "User-Agent" => "Mozilla/5.0 (X11; Linux x86_64; rv:120.0) Gecko/20100101 Firefox/120.0",
        # "Content-Length" => "229", # that should not be fixed
    )
    # BODY, that's the part that changes
    # format the time
    datetime = DateTime(year, month, day, hour, minute)
    datetime = Dates.format(datetime, "yyyy-mm-ddTHH:MM:SS.000Z")
    body = Dict(
        "timeType" => "ut",         # universal time
        "coordinateType" => "geo",  # geographic
        "datetime" => datetime,
        "lat" => lat,
        "lon" => lon,
        "height" => 100,
        "profileType" => 1,       # height profile
        "start" => height[1],
        "stop" => height[end],
        "step" => step(height),

        # Everything under are default parameters
        "tecUpper" => 2000,
        "tecLower" => 50,
        "outputType" => 0,
        "useOptionals" => false,
        "layVersion" => false,
        "NeTopside" => "NeQuick",
        "fof2Model" => "URSI-88",
        "fof2Storm" => true,
        "NeTopsideStorm" => false,
        "hmF2Model"=> "AMTB-model",
        "bottomsideThicknessB0" => "ABT-2009",
        "F1Model" => "Scotto-1997-no-L",
        "EPeakAuroralStorm" => false,
        "D" => "IRI-1990",
        "Te" => "TBT-2012_PF107",
        "Ti" => "Tru-2021",
        "ionCompModel" => "RBV10/TBT15",
        "auroralBoundaryModel" => false,
    )
    body = JSON.json(body) # convert to JSON
    output_file_link = HTTP.request("POST", url; headers, body, connect_timeout = 300)     # request results
    output_file_link = JSON.parse(String(output_file_link.body))    # from JSON to Julia
    output_file_link = output_file_link["txt"]                      # get only the link


    # read the results
    iri_data = readdlm(Downloads.download(output_file_link))
    println("Iri data succesfully downloaded.")


    parameters = (; year, month, day, hour, minute, lat, lon, height)
    return iri_data, parameters
end



using PythonCall
function calculate_iri_data(year=2018, month=12, day=7, hour=11, minute=15,
    lat = 76, lon = 5, height = 85:1:700)

    print("Calculating iri data...")

    datetime = pyimport("datetime")
    time = datetime.datetime(year, month, day, hour, minute, 0)

    # import iri2016 model from the Python package 'iri2016'
    iri2016 = pyimport("iri2016.profile")
    # run the model
    iri_data = iri2016.IRI(time, Py([height[1], height[end], step(height)]), Py(lat), Py(lon))
    # convert the Python Dataset to a DataArray
    iri_data = iri_data.to_dataarray()
    # convert from Python array to Julia array
    iri_data = pyconvert(Array, iri_data) # of size (20, n_z, 1)
    # change shape from (20, n_z, 1) to (n_z, 20)
    iri_data = iri_data[:, :, 1]'
    # add a column with the altitude
    iri_data = hcat(Vector(height), iri_data)
    # add a header with the name of columns
    iri_data = vcat(["height(km)" "ne(m-3)" "Tn(K)" "Ti(K)" "Te(K)" "nO+(m-3)" "nH+(m-3)" "nHe+(m-3)" "nO2+(m-3)" "nNO+(m-3)" "nCI(m-3)" "nN+(m-3)" "NmF2" "hmF2" "NmF1" "hmF1" "NmE" "hmE" "TEC" "EqVertIonDrift" "foF2"], iri_data)

    println(" done.")

    parameters = (; year, month, day, hour, minute, lat, lon, height)
    return iri_data, parameters
end


function save_iri_data(iri_data, parameters)
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
        writedlm(f, iri_data)
    end
    println("File " * @bold("$filename") * " created under " * @underline("$directory") * ".")

    return fullpath
end

function load_parameters_iri(iri_file)
    local data
    # Loading the `parameters` section of the file
    open(`head -n12 $iri_file`) do io
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
