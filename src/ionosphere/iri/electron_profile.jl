using Dates: DateTime

# ======================================================================================== #
#                           ElectronProfile (electron source)                              #
# ======================================================================================== #

"""
    ElectronProfile(h, ne, Te; origin="")
    ElectronProfile{T}(h, ne, Te; origin="")

Electron density `ne` (m⁻³) and temperature `Te` (K) on the altitude grid `h` (m). Callable on
any altitude grid (m); returns `(; ne, Te)` interpolated to that grid (`ne` in log-space, `Te`
linearly). Produced by [`run_iri`](@ref), [`read_iri_file`](@ref) and [`read_ccmc_iri`](@ref).

`T` is the floated promotion of the input element types; `ElectronProfile{T}` converts the
inputs to `T` instead. `origin` is a free-form provenance label, shown by `show` and written
into `inputs/atmosphere.nc`.

# Example
```julia
profile = ElectronProfile(h_m, ne_m3, Te_K; origin="my measurement")
ne, Te  = profile(altitude_grid.h)
```
"""
struct ElectronProfile{T<:Real}
    h::Vector{T}          # native altitude (m)
    ne::Vector{T}         # electron density (m⁻³)
    Te::Vector{T}         # electron temperature (K)
    origin::String        # provenance label (free-form, may be empty)

    function ElectronProfile{T}(h, ne, Te, origin) where {T<:Real}
        # Convert before validating so the checks see the values that will be stored.
        h  = convert(Vector{T}, h)
        ne = convert(Vector{T}, ne)
        Te = convert(Vector{T}, Te)
        check_profile_grid("ElectronProfile", h, ("ne", ne), ("Te", Te))
        return new{T}(h, ne, Te, String(origin))
    end
end

ElectronProfile{T}(h, ne, Te; origin::AbstractString = "") where {T<:Real} =
    ElectronProfile{T}(h, ne, Te, origin)

ElectronProfile(h, ne, Te, origin) =
    ElectronProfile{promote_type(float(eltype(h)), float(eltype(ne)),
                                 float(eltype(Te)))}(h, ne, Te, origin)
ElectronProfile(h, ne, Te; origin::AbstractString = "") =
    ElectronProfile(h, ne, Te, origin)

function (p::ElectronProfile)(h_atm::AbstractVector)
    warn_extrapolation(p, h_atm)
    ne = interpolate_profile(p.ne, p.h ./ 1e3, h_atm; log_interpolation = true)
    Te = interpolate_profile(p.Te, p.h ./ 1e3, h_atm; log_interpolation = false)
    return (; ne, Te)
end

Base.show(io::IO, p::ElectronProfile) = print(io, profile_label(p))

function Base.show(io::IO, ::MIME"text/plain", p::ElectronProfile)
    println(io, "ElectronProfile:")
    println(io, "├── Origin:    ", isempty(p.origin) ? "(unlabelled)" : p.origin)
    println(io, "├── Altitudes: ", length(p.h),
                " ($(p.h[1] / 1e3) – $(p.h[end] / 1e3) km)")
    println(io, "├── Max ne:    ", round(maximum(p.ne), sigdigits=3), " m⁻³")
    print(io,   "└── Max Te:    ", round(maximum(p.Te), sigdigits=3), " K")
end


# ======================================================================================== #
#                              Producers                                                   #
# ======================================================================================== #

"""
    run_iri(; year=2018, month=12, day=7, hour=11, minute=15, lat=76, lon=5,
             height=85:1:700, save_to=nothing, verbose=true) -> ElectronProfile

Run the IRI-2020 model (Python `iri2020` package) for the given conditions and return `ne` and
`Te` as an [`ElectronProfile`](@ref). Levels where IRI reports the -1 sentinel (no valid
profile, typically at the bottom of the range) are dropped with a warning.

# Keyword Arguments
- `height`: altitude levels (km) at which the model is evaluated.
- `save_to`: directory in which to also write the model output as an AURORA IRI text file,
  readable with [`read_iri_file`](@ref). Saving into `internal_data/data_electron/` makes the
  file visible to [`find_iri_file`](@ref).

# Example
```julia
electrons = run_iri(; year=2005, month=10, day=8, hour=22, minute=0, lat=69.58, lon=19.23)
model = AuroraModel(altitude_lims, θ_lims, E_max, neutrals, electrons)
```
"""
function run_iri(; year = 2018, month = 12, day = 7, hour = 11, minute = 15,
                 lat = 76, lon = 5, height = 85:1:700,
                 save_to = nothing, verbose = true)
    iri_data, parameters = calculate_iri_data(; year, month, day, hour, minute, lat, lon,
                                               height, verbose)
    if save_to !== nothing
        save_iri_data(iri_data, parameters; directory = save_to, verbose)
    end
    data        = iri_data[2:end, :]               # drop the header row
    instant     = DateTime(year, month, day, hour, minute)
    origin      = "IRI-2020 $instant $(lat)N/$(lon)E"

    raw = (height_km = Float64.(data[:, 1]),
           ne        = Float64.(data[:, 2]),       # electron density (m⁻³)
           Te        = Float64.(data[:, 5]))       # electron temperature (K)
    trimmed = trim_iri_sentinels(raw, "IRI-2020 run for $instant at $(lat)N/$(lon)E\n")

    return ElectronProfile(trimmed.height_km .* 1e3, trimmed.ne, trimmed.Te; origin)
end

"""
    read_iri_file(iri_file) -> ElectronProfile

Read `ne` and `Te` from an IRI text file generated by AURORA (see [`find_iri_file`](@ref))
and return them as an [`ElectronProfile`](@ref) on the file's altitude grid.
"""
function read_iri_file(iri_file::AbstractString)
    raw = load_iri(iri_file)
    return ElectronProfile(raw.data.height_km .* 1e3, raw.data.ne, raw.data.Te;
                           origin = "IRI file $(basename(iri_file))")
end

"""
    read_ccmc_iri(file) -> ElectronProfile

Read `ne` and `Te` from a CCMC ModelWeb IRI text export and return them as an
[`ElectronProfile`](@ref) (converted to m⁻³, `-1` sentinel levels dropped). The columns
`km`, `Ne/cm-3` and `Te/K` are located by header name; a missing one is an error.

# Example
```julia
electrons = read_ccmc_iri("iri_output.txt")
model = AuroraModel(altitude_lims, θ_lims, E_max, neutrals, electrons)
```
"""
function read_ccmc_iri(file::AbstractString)
    lines = readlines(file)
    # Detect the header on the unit-free "Ne/" (plus a "km" token, so preamble prose does not
    # match), then check the exact names. An export in other units fails here.
    header_idx, header, column, columns_found = locate_ccmc_header(
        lines, l -> occursin("Ne/", l) && "km" in split(l), file, "read_ccmc_iri",
        "a line naming an \"Ne/…\" column")
    for name in ("km", "Ne/cm-3", "Te/K")
        haskey(column, name) || throw(ArgumentError(
            "read_ccmc_iri: no \"$name\" column in the header of $file.\n" * columns_found))
    end
    h_col, ne_col, Te_col = column["km"], column["Ne/cm-3"], column["Te/K"]
    n_cols = max(h_col, ne_col, Te_col)

    h_km = Float64[]
    ne   = Float64[]
    Te   = Float64[]
    for l in lines[(header_idx + 1):end]
        cols = split(l)
        length(cols) >= n_cols || continue
        h  = tryparse(Float64, cols[h_col])
        n  = tryparse(Float64, cols[ne_col])
        t  = tryparse(Float64, cols[Te_col])
        (h === nothing || n === nothing || t === nothing) && continue
        # Drop unusable levels. CCMC marks them with -1, but ne is interpolated in log-space
        # downstream, so zero is just as fatal and is dropped the same way.
        (n <= 0 || t <= 0) && continue
        push!(h_km, h); push!(ne, n); push!(Te, t)
    end
    isempty(h_km) && throw(ArgumentError(
        "read_ccmc_iri: no valid (non-sentinel) data rows parsed from $file"))

    return ElectronProfile(h_km .* 1e3, ne .* 1e6, Te;   # km→m, cm⁻³→m⁻³
                           origin = "CCMC IRI $(basename(file))")
end
