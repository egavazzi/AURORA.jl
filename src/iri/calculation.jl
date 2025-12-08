using PythonCall: pyimport, pyconvert, Py

"""
    calculate_iri_data(; year=2018, month=12, day=7, hour=11, minute=15,
                        lat=76, lon=5, height=85:1:700)

Calculate IRI-2016 ionospheric model data using Python interface.

This function calls the Python `iri2016` package to compute ionospheric parameters
including electron density, temperatures, and ion composition profiles. The data
is returned as a matrix with a header row containing column names.

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
  - Matrix with header row and IRI data (21 columns: height + 20 parameters)
  - NamedTuple with input parameters

# Data Columns
The returned matrix contains the following columns:
1. height(km): Altitude
2. ne(m⁻³): Electron density
3. Tn(K): Neutral temperature
4. Ti(K): Ion temperature
5. Te(K): Electron temperature
6-12. nO⁺, nH⁺, nHe⁺, nO2⁺, nNO⁺, nCI, nN⁺: Ion densities (m⁻³)
13-18. NmF2, hmF2, NmF1, hmF1, NmE, hmE: Peak densities and heights
19. TEC: Total Electron Content
20. EqVertIonDrift: Equatorial vertical ion drift
21. foF2: F2 critical frequency
"""
function calculate_iri_data(; year = 2018, month = 12, day = 7, hour = 11, minute = 15,
                             lat = 76, lon = 5, height = 85:1:700)
    print("Calculating iri data...")
    datetime = pyimport("datetime")
    time = datetime.datetime(year, month, day, hour, minute, 0)
    # import iri2016 model from the Python package 'iri2016'
    redirect_stdio(stdout=devnull) do # hide all the extra verbosity from iri2016
        global iri2016 = pyimport("iri2016.profile")
    end
    # iri2016 = pyimport("iri2016.profile")
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
    parameters = (; year, month, day, hour, minute, lat, lon, height)
    println(" done.")
    return iri_data, parameters
end
