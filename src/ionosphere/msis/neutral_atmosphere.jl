using Dates: DateTime

# ======================================================================================== #
#                              Density-source types                                        #
# ======================================================================================== #

"""
    DensityProfile(h, n; origin="")
    DensityProfile{T}(h, n; origin="")

Density `n` (m⁻³) on the altitude grid `h` (m). Callable on any altitude grid (m); evaluates
via PCHIP interpolation in log-space. This is the density source of a [`NeutralSpecies`](@ref)
and the per-species content of a [`NeutralAtmosphere`](@ref).

`T` is the floated promotion of the input element types; `DensityProfile{T}` converts the
inputs to `T` instead. `origin` is a free-form provenance label, shown by `show` and written
into `inputs/atmosphere.nc`.

# Example
```julia
profile = DensityProfile(h_msis_m, n_N2; origin="ccmc_run_4321.txt")
n = profile(altitude_grid.h)
```
"""
struct DensityProfile{T<:Real}
    h::Vector{T}         # altitude (m)
    n::Vector{T}         # density (m⁻³)
    origin::String       # provenance label (free-form, may be empty)

    function DensityProfile{T}(h, n, origin) where {T<:Real}
        # Convert before validating so the checks see the values that will be stored.
        h = convert(Vector{T}, h)
        n = convert(Vector{T}, n)
        check_profile_grid("DensityProfile", h, ("n", n))
        return new{T}(h, n, String(origin))
    end
end

DensityProfile{T}(h, n; origin::AbstractString = "") where {T<:Real} =
    DensityProfile{T}(h, n, origin)

DensityProfile(h, n, origin) =
    DensityProfile{promote_type(float(eltype(h)), float(eltype(n)))}(h, n, origin)
DensityProfile(h, n; origin::AbstractString = "") = DensityProfile(h, n, origin)

function (d::DensityProfile)(h_atm::AbstractVector)
    warn_extrapolation(d, h_atm)
    return interpolate_profile(d.n, d.h ./ 1e3, h_atm; log_interpolation = true)
end

Base.show(io::IO, d::DensityProfile) = print(io, profile_label(d))

function Base.show(io::IO, ::MIME"text/plain", d::DensityProfile)
    println(io, "DensityProfile:")
    println(io, "├── Origin:    ", isempty(d.origin) ? "(unlabelled)" : d.origin)
    println(io, "├── Altitudes: ", length(d.h),
                " ($(d.h[1] / 1e3) – $(d.h[end] / 1e3) km)")
    print(io,   "└── Max n:     ", round(maximum(d.n), sigdigits=3), " m⁻³")
end


# ======================================================================================== #
#                              NeutralAtmosphere                                           #
# ======================================================================================== #

"""
    NeutralAtmosphere(densities; origin="", dropped=Symbol[])
    NeutralAtmosphere(:N2 => profile, :O2 => profile, ...; origin="", dropped=Symbol[])

One [`DensityProfile`](@ref) per species, keyed by symbol (`:N2`, `:O2`, `:O`, `:He`, `:H`,
`:Ar`, `:N`, `:NO`). Pass it as the `neutrals` argument of [`AuroraModel`](@ref), or index it
(`neutrals[:N2]`) to get one species' density source. Produced by [`run_msis`](@ref),
[`read_msis_file`](@ref) and [`read_ccmc_msis`](@ref), or built from a dictionary or from
`species => profile` pairs.

Read-only collection interface: `getindex`, `haskey`, `get`, `keys`, `values`, `pairs`,
`length`, and iteration over `species => profile` pairs. Each species keeps its own altitude
grid, restricted to the levels where the source reports it.

# Example
```julia
neutrals = read_ccmc_msis("nrlmsis_output.txt")
model    = AuroraModel(altitude_lims, θ_lims, E_max, neutrals, electrons)
n_N2     = neutrals[:N2](altitude_grid.h)

neutrals = NeutralAtmosphere(:N2 => DensityProfile(h, n_N2), :O2 => DensityProfile(h, n_O2);
                             origin = "my radar inversion")
```
"""
struct NeutralAtmosphere
    # Abstract value type on purpose: species may carry different element types.
    densities::Dict{Symbol, DensityProfile}
    origin::String
    dropped::Vector{Symbol}   # species the source carried but never usably reported

    function NeutralAtmosphere(densities, origin, dropped)
        return new(Dict{Symbol, DensityProfile}(densities), String(origin),
                   collect(Symbol, dropped))
    end
end

NeutralAtmosphere(densities::AbstractDict; origin::AbstractString = "",
                  dropped = Symbol[]) =
    NeutralAtmosphere(densities, origin, dropped)

NeutralAtmosphere(first_pair::Pair{Symbol, <:DensityProfile},
                  rest::Pair{Symbol, <:DensityProfile}...;
                  origin::AbstractString = "", dropped = Symbol[]) =
    NeutralAtmosphere(Dict{Symbol, DensityProfile}(first_pair, rest...), origin, dropped)

function Base.getindex(p::NeutralAtmosphere, species::Symbol)
    if !haskey(p.densities, species)
        hint = species in p.dropped ?
            " It is present in the source but reported at fewer than 2 usable levels." : ""
        throw(ArgumentError(
            "NeutralAtmosphere: no density for :$species.$hint Available: " *
            join(sort!(string.(keys(p.densities))), ", ")))
    end
    return p.densities[species]
end

Base.haskey(p::NeutralAtmosphere, species::Symbol) = haskey(p.densities, species)
Base.get(p::NeutralAtmosphere, species::Symbol, default) = get(p.densities, species, default)
Base.keys(p::NeutralAtmosphere) = keys(p.densities)
Base.values(p::NeutralAtmosphere) = values(p.densities)
Base.pairs(p::NeutralAtmosphere) = pairs(p.densities)
Base.length(p::NeutralAtmosphere) = length(p.densities)
Base.eltype(::Type{NeutralAtmosphere}) = Pair{Symbol, DensityProfile}
Base.iterate(p::NeutralAtmosphere) = iterate(p.densities)
Base.iterate(p::NeutralAtmosphere, state) = iterate(p.densities, state)

function Base.show(io::IO, p::NeutralAtmosphere)
    print(io, "NeutralAtmosphere(", join(sort!(string.(keys(p.densities))), ", "), ")")
end

function Base.show(io::IO, ::MIME"text/plain", p::NeutralAtmosphere)
    println(io, "NeutralAtmosphere:")
    println(io, "├── Origin:  ", isempty(p.origin) ? "(unlabelled)" : p.origin)
    print(io,   "└── Species: ", join(sort!(string.(keys(p.densities))), ", "))
end

# Number of leading header names that line up with their data column. CCMC's NRLMSISE-00
# export writes "Heden(cm-3)Arden(cm-3)" without a space, so names from that token on would
# address the wrong data column. The data column count is the modal token count over the
# numeric rows, so that one malformed row does not decide the shape.
function trusted_header_width(header, lines, header_idx, file)
    token_counts = Int[]
    for l in @view lines[(header_idx + 1):end]
        isempty(strip(l)) && continue
        tokens = split(l)
        tryparse(Float64, tokens[1]) === nothing && continue
        push!(token_counts, length(tokens))
    end
    isempty(token_counts) && return length(header)
    n_data = argmax(n -> count(==(n), token_counts), unique(token_counts))
    n_data == length(header) && return length(header)

    # A name with ')' before its last character is two names glued together.
    glued = findfirst(t -> occursin(r"\)\S", t), header)
    width = glued === nothing ? 0 : glued - 1
    loc   = glued === nothing ? "" :
            ", starting at the run-together name '" * header[glued] * "'"
    @warn "read_ccmc_msis: $(basename(file)) has $(length(header)) header names for " *
          "$n_data data columns$loc. Columns from there on cannot be matched to their " *
          "data and are ignored; the species named before it are read normally." maxlog = 1
    return width
end

# Levels where a species is reported. Unreported levels arrive as NaN (pymsis writes NaN;
# other readers translate their own marker); non-positive values cannot be log-interpolated.
usable_levels(n) = findall(x -> isfinite(x) && x > 0, n)

# Build the per-species density sources of a NeutralAtmosphere, keeping for each species only
# the levels where it is reported. `densities` maps a species to its column on the grid `h_m`.
function species_densities(h_m, densities, label)
    out     = Dict{Symbol, DensityProfile}()
    dropped = Symbol[]
    for (species, n) in densities
        valid = usable_levels(n)
        if length(valid) < 2
            push!(dropped, species)
            continue
        end
        out[species] = DensityProfile(h_m[valid], n[valid]; origin = "$label :$species")
    end
    return out, dropped
end


# ======================================================================================== #
#                              Producers                                                   #
# ======================================================================================== #

"""
    run_msis(; year=2018, month=12, day=7, hour=11, minute=15, lat=76, lon=5,
              height=85:1:700, save_to=nothing, verbose=true) -> NeutralAtmosphere

Run the NRLMSIS 2.1 model (Python `pymsis` package) for the given conditions and return the
species densities as a [`NeutralAtmosphere`](@ref). Levels where the model does not report a
species (`NaN`, e.g. N at low altitude) are dropped for that species only.

# Keyword Arguments
- `height`: altitude levels (km) at which the model is evaluated.
- `save_to`: directory in which to also write the model output as an AURORA MSIS text file,
  readable with [`read_msis_file`](@ref). Saving into `internal_data/data_neutrals/` makes the
  file visible to [`find_msis_file`](@ref).

# Example
```julia
neutrals  = run_msis(; year=2005, month=10, day=8, hour=22, minute=0, lat=69.58, lon=19.23)
electrons = run_iri(; year=2005, month=10, day=8, hour=22, minute=0, lat=69.58, lon=19.23)
model     = AuroraModel(altitude_lims, θ_lims, E_max, neutrals, electrons)
```
"""
function run_msis(; year = 2018, month = 12, day = 7, hour = 11, minute = 15,
                  lat = 76, lon = 5, height = 85:1:700,
                  save_to = nothing, verbose = true)
    msis_data, parameters = calculate_msis_data(; year, month, day, hour, minute, lat, lon,
                                                 height, verbose)
    if save_to !== nothing
        save_msis_data(msis_data, parameters; directory = save_to, verbose)
    end
    data    = msis_data[2:end, :]                  # drop the header row
    h_m     = Float64.(data[:, 1]) .* 1e3          # height(km) → m
    instant = DateTime(year, month, day, hour, minute)
    label   = "NRLMSIS 2.1 $instant $(lat)N/$(lon)E"

    # Columns of calculate_msis_data, all number densities in m⁻³
    column_spec = (:N2 => 3, :O2 => 4, :O => 5, :He => 6, :H => 7, :Ar => 8, :N => 9,
                   :NO => 11)
    columns     = Dict(s => Float64.(data[:, c]) for (s, c) in column_spec)

    densities, dropped = species_densities(h_m, columns, label)
    return NeutralAtmosphere(densities; origin = label, dropped)
end

"""
    read_msis_file(msis_file) -> NeutralAtmosphere

Read every species from an MSIS text file generated by AURORA (see [`find_msis_file`](@ref))
and return them as a [`NeutralAtmosphere`](@ref) on the file's altitude grid.
"""
function read_msis_file(msis_file::AbstractString)
    raw   = load_msis(msis_file)
    h_m   = raw.data.height_km .* 1e3
    label = "MSIS file $(basename(msis_file))"

    columns = Dict(species => getproperty(raw.data, species)
                   for species in (:N2, :O2, :O, :He, :H, :Ar, :N, :NO)
                   if hasproperty(raw.data, species))
    densities, dropped = species_densities(h_m, columns, label)
    return NeutralAtmosphere(densities; origin = label, dropped)
end

"""
    read_ccmc_msis(file) -> NeutralAtmosphere

Read the species densities from a CCMC ModelWeb NRLMSIS text export (NRLMSIS 2.x or
NRLMSISE-00) and return them as a [`NeutralAtmosphere`](@ref), converted to m⁻³, with the
`9.999E-38` sentinel levels dropped per species. Columns (`Heit(km)`, `N2den(cm-3)`, …) are
located by header name; species absent from the export are not returned.

# Example
```julia
neutrals = read_ccmc_msis("nrlmsis_output.txt")
model    = AuroraModel(altitude_lims, θ_lims, E_max, neutrals, electrons)
```
"""
function read_ccmc_msis(file::AbstractString)
    lines = readlines(file)
    # Columns are matched by full header name, unit included, so a change of unit is an error.
    header_idx, header, column, columns_found = locate_ccmc_header(
        lines, l -> occursin("N2den", l), file, "read_ccmc_msis", "a line containing \"N2den\"")
    trusted = trusted_header_width(header, lines, header_idx, file)

    haskey(column, "Heit(km)") || throw(ArgumentError(
        "read_ccmc_msis: no altitude column \"Heit(km)\" in the header of $file.\n" *
        columns_found))
    h_col = column["Heit(km)"]
    h_col <= trusted || throw(ArgumentError(
        "read_ccmc_msis: the altitude column \"Heit(km)\" of $file sits past a run-together " *
        "header name, so it cannot be matched to its data column.\n" * columns_found))

    species_columns = [s => column[name] for (s, name) in
                       (:O  => "Oden(cm-3)",  :N2 => "N2den(cm-3)", :O2 => "O2den(cm-3)",
                        :NO => "NOden(cm-3)", :He => "Heden(cm-3)", :Ar => "Arden(cm-3)",
                        :H  => "Hden(cm-3)",  :N  => "Nden(cm-3)")
                       if haskey(column, name) && column[name] <= trusted]
    isempty(species_columns) && throw(ArgumentError(
        "read_ccmc_msis: no known species density column (\"N2den(cm-3)\", \"Oden(cm-3)\", " *
        "…) in the header of $file.\n" * columns_found))

    n_cols = max(h_col, maximum(last, species_columns))
    h_km   = Float64[]
    raw    = Dict(s => Float64[] for (s, _) in species_columns)
    for l in lines[(header_idx + 1):end]
        cols = split(l)
        length(cols) >= n_cols || continue
        h = tryparse(Float64, cols[h_col])
        h === nothing && continue
        # An unparseable field drops this level for that species only.
        push!(h_km, h)
        for (s, c) in species_columns
            v = tryparse(Float64, cols[c])
            push!(raw[s], v === nothing ? NaN : v)
        end
    end
    isempty(h_km) && throw(ArgumentError(
        "read_ccmc_msis: no valid data rows parsed from $file"))

    label = "CCMC NRLMSIS $(basename(file))"
    # CCMC writes 9.999E-38 where it does not report a species. Turn those into NaN so they
    # are dropped per species, and convert the rest cm⁻³ → m⁻³.
    columns = Dict(species => [v > 1e-37 ? v * 1e6 : NaN for v in raw[species]]
                   for (species, _) in species_columns)
    densities, dropped = species_densities(h_km .* 1e3, columns, label)
    return NeutralAtmosphere(densities; origin = label, dropped)
end
