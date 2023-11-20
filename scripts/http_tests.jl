using HTTP, JSON

url = "https://kauai.ccmc.gsfc.nasa.gov/instantrun/api/nrlmsis"

r = HTTP.request("GET", "https://kauai.ccmc.gsfc.nasa.gov/instantrun/api/nrlmsis")

headers = Dict(
    "Host" => "kauai.ccmc.gsfc.nasa.gov",
    "User-Agent" => "Mozilla/5.0 (X11; Linux x86_64; rv:109.0) Gecko/20100101 Firefox/119.0",
    "Accept" => "application/json",
    "Accept-Language" => "no,sv;q=0.8,fr;q=0.6,en-GB;q=0.4,en;q=0.2",
    "Accept-Encoding" => "gzip, deflate, br",
    "Referer" => "https://kauai.ccmc.gsfc.nasa.gov/instantrun/nrlmsis/",
    "Content-Length" => "229",
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

body = "{\"version\":\"NRLMSIS2\",\"timeType\":\"ut\",\"coordinateType\":\"geo\",\"datetime\":\"2023-11-14T12:00:00.000Z\",\"lat\":80,\"lon\":45,\"height\":100,\"profileType\":1,\"start\":0,\"stop\":1000,\"step\":50,\"fdaily\":-1,\"f3month\":-1,\"ap\":-1,\"_datelimit\":true}"
JSON.parse(body)
# body = JSON.json(body)




output_file_link = HTTP.request("POST", url, headers, body) # request results
output_file_link = JSON.parse(String(output_file_link.body))    # from JSON to Julia
output_file_link = output_file_link["txt"]

using DelimitedFiles
using Downloads

readdlm(Downloads.download(output_file_link), skipstart = 1)
