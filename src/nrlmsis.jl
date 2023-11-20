using HTTPS, JSON, Downloads, DelimitedFiles

function download_msis_data(;
        year = 2018,
        month = 12,
        day = 7,
        hour = 11,
        minute = 15,
        lat = 76,
        lon = 5,
        height = 85:1:700,
    )
    # Initialize natrix
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
        println(local_height)
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
        output_file_link = HTTP.request("POST", url, headers, body)     # request results
        output_file_link = JSON.parse(String(output_file_link.body))    # from JSON to Julia
        output_file_link = output_file_link["txt"]                      # get only the link

        if i_query == 1
            # first one, we keep the headers (skipstart = 1)
            nrlmsis_data = readdlm(Downloads.download(output_file_link), skipstart = 1)
        else
            # we skip the headers (skipstart = 2)
            data = readdlm(Downloads.download(output_file_link), skipstart = 2)
            nrlmsis_data = vcat(nrlmsis_data, data)
        end
    end

    return nrlmsis_data
end

function save_msis_data(year, month, day, hour, minute, lat, lon, height, nrlmsis_data)

    filename = "msis_$year$month$day-$(hour)$(minute)_$(lat)N$(lon)E.txt"
    # println(filename)

end
    directory = pkgdir(AURORA, "internal_data", "data_neutrals")



save_msis_data(2022, 10, 10, 12, 10, 70, 80, 80:1:500, 10)



download_msis_data(height = 85:1:700)




# a = 2018
# b = 12
# c = 7
# d = 11
# e = 15
# lat = 76
# lon = 5
# height = 85:1:700
