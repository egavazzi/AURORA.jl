using PythonCall: pyimport, pyconvert, Py
using Downloads: Downloads
using DelimitedFiles: writedlm

"""
    calculate_msis_data(; year=2018, month=12, day=7, hour=11, minute=15,
                         lat=76, lon=5, height=85:1:700)

Calculate NRLMSIS 2.1 atmospheric model data using Python interface.

This function calls the Python `pymsis` package to compute neutral atmosphere
parameters including densities of various species (N₂, O₂, O, He, H, Ar, N, NO)
and temperature profiles. The data is returned as a matrix with a header row
containing column names.

# Keyword Arguments
- `year::Int=2018`: Year (defaults to Visions2 launch conditions)
- `month::Int=12`: Month (1-12)
- `day::Int=7`: Day of month (1-31)
- `hour::Int=11`: Hour in Universal Time (0-23)
- `minute::Int=15`: Minute (0-59)
- `lat::Real=76`: Geographic latitude in degrees North
- `lon::Real=5`: Geographic longitude in degrees East
- `height::AbstractRange=85:1:700`: Altitude range in km

# Returns
- `Tuple{Matrix, NamedTuple}`:
  - Matrix with header row and MSIS data (12 columns: height + 11 parameters)
  - NamedTuple with input parameters

# Data Columns
The returned matrix contains the following columns:
1. height(km): Altitude
2. air(kg/m³): Total mass density
3. N₂(m⁻³): Molecular nitrogen density
4. O₂(m⁻³): Molecular oxygen density
5. O(m⁻³): Atomic oxygen density
6. He(m⁻³): Helium density
7. H(m⁻³): Atomic hydrogen density
8. Ar(m⁻³): Argon density
9. N(m⁻³): Atomic nitrogen density
10. anomalousO(m⁻³): Anomalous oxygen density
11. NO(m⁻³): Nitric oxide density
12. T(K): Temperature

# Notes
- This function automatically downloads the SW-All.csv file (space weather data)
  if not already present in the pymsis package directory
- The geomagnetic activity parameter is set to -1 (use observational data)
"""
function calculate_msis_data(; year = 2018, month = 12, day = 7, hour = 11, minute = 15,
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
    SW_file = joinpath(path_to_pymsis, "SW-All.csv")
    if isfile(SW_file) == false
        print(" downloading ap and F10.7 data from (https://celestrak.org/SpaceData/SW-All.csv)...")
        Downloads.download("https://celestrak.org/SpaceData/SW-All.csv", SW_file)
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
