using DataInterpolations: PCHIPInterpolation, ExtrapolationType
using Dates: DateTime

# ======================================================================================== #
#                              Density-source types                                        #
# ======================================================================================== #

"""
    DensityProfile(h, n; source="")

Density source defined by user-supplied altitude (`h`, m) and density (`n`, m⁻³) vectors.
Callable on any altitude grid (m); evaluates via PCHIP interpolation in log-space,
consistent with AURORA's MSIS interpolation convention.

This is the universal interchange for densities produced outside AURORA (CCMC ModelWeb runs,
radar inversions, any external atmospheric model): reduce the source to an altitude vector and
a density vector, then wrap it here. The optional `source` string records provenance; it is
shown by `show` and written into `inputs/atmosphere.nc`.

# Example
```julia
profile = DensityProfile(h_msis_m, n_N2; source="ccmc_run_4321.txt")
n = profile(altitude_grid.h)
```
"""
struct DensityProfile
    h::Vector{Float64}   # altitude (m)
    n::Vector{Float64}   # density (m⁻³)
    source::String       # provenance label (free-form, may be empty)
end

function DensityProfile(h, n; source::AbstractString = "")
    h, n = collect(Float64, h), collect(Float64, n)
    check_profile_grid("DensityProfile", h, ("n", n))
    return DensityProfile(h, n, String(source))
end

function (d::DensityProfile)(h_atm::AbstractVector)
    itp = PCHIPInterpolation(log.(d.n), d.h ./ 1e3;
                             extrapolation = ExtrapolationType.Linear)
    return exp.(itp(collect(Float64, h_atm) ./ 1e3))
end

Base.show(io::IO, d::DensityProfile) = print(io, profile_label(d))

function Base.show(io::IO, ::MIME"text/plain", d::DensityProfile)
    println(io, "DensityProfile:")
    println(io, "├── Source:    ", isempty(d.source) ? "(unlabelled)" : d.source)
    println(io, "├── Altitudes: ", length(d.h),
                " ($(d.h[1] / 1e3) – $(d.h[end] / 1e3) km)")
    print(io,   "└── Max n:     ", round(maximum(d.n), sigdigits=3), " m⁻³")
end


# ======================================================================================== #
#                           NeutralProfile (neutral atmosphere)                             #
# ======================================================================================== #

"""
    NeutralProfile(densities; source="")

Neutral atmosphere holding one [`DensityProfile`](@ref) per species, keyed by symbol
(`:N2`, `:O2`, `:O`, `:He`, `:H`, `:Ar`, `:N`, `:NO`). Index it to get a single species'
density source (`neutrals[:N2]`), or pass it directly as the `neutrals` argument of
[`AuroraModel`](@ref) to build the three default species from it.

This is the neutral analogue of [`ElectronProfile`](@ref): the universal interchange for a
full set of neutral densities, whatever their origin. Build one from a legacy AURORA MSIS
file with [`read_msis_file`](@ref), from a CCMC ModelWeb NRLMSIS download with
[`read_ccmc_msis`](@ref), or directly from your own vectors. Because it stores data (not a
file path), it round-trips through `physics_state.jld2` and reproduces on any machine with
no external file.

Species are stored on their own native altitude grids, so a source that reports a species
only over part of the column keeps just the levels where that species is actually defined.
Each reader translates its own missing-value marker before building the profile.

# Example
```julia
neutrals = read_ccmc_msis("nrlmsis_output.txt")
model    = AuroraModel(altitude_lims, θ_lims, E_max, neutrals, electrons)

n_N2 = neutrals[:N2](altitude_grid.h)   # sample one species directly
```
"""
struct NeutralProfile
    densities::Dict{Symbol, DensityProfile}
    source::String
    dropped::Vector{Symbol}   # species the source carried but never usably reported
end

NeutralProfile(densities::AbstractDict; source::AbstractString = "",
               dropped = Symbol[]) =
    NeutralProfile(Dict{Symbol, DensityProfile}(densities), String(source),
                   collect(Symbol, dropped))

function Base.getindex(p::NeutralProfile, species::Symbol)
    if !haskey(p.densities, species)
        # Distinguish "your source never mentioned this" from "it did, but reported nothing
        # usable" — otherwise a species silently dropped at read time looks like a typo.
        hint = species in p.dropped ?
            " It is present in the source but reported at fewer than 2 usable levels." : ""
        throw(ArgumentError(
            "NeutralProfile: no density for :$species.$hint Available: " *
            join(sort!(string.(keys(p.densities))), ", ")))
    end
    return p.densities[species]
end

Base.haskey(p::NeutralProfile, species::Symbol) = haskey(p.densities, species)
Base.keys(p::NeutralProfile) = keys(p.densities)

function Base.show(io::IO, p::NeutralProfile)
    print(io, "NeutralProfile(", join(sort!(string.(keys(p.densities))), ", "), ")")
end

function Base.show(io::IO, ::MIME"text/plain", p::NeutralProfile)
    println(io, "NeutralProfile:")
    println(io, "├── Source:  ", isempty(p.source) ? "(unlabelled)" : p.source)
    print(io,   "└── Species: ", join(sort!(string.(keys(p.densities))), ", "))
end

# Indices of the levels where a species is actually reported. Anything a source could not
# give us is expected to arrive as NaN — pymsis (and so AURORA's own MSIS files) writes NaN
# directly, while a format with its own missing-value marker translates it before getting
# here. None of these can go through the log-space interpolation.
usable_levels(n) = findall(x -> isfinite(x) && x > 0, n)

# Build the per-species density sources of a NeutralProfile, keeping for each species only
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
        out[species] = DensityProfile(h_m[valid], n[valid]; source = "$label :$species")
    end
    return out, dropped
end


# ======================================================================================== #
#                              Producers                                                   #
# ======================================================================================== #

"""
    run_msis(; year=2018, month=12, day=7, hour=11, minute=15, lat=76, lon=5,
              height=85:1:700, verbose=true) -> NeutralProfile

Run the NRLMSIS 2.1 model (via the Python `pymsis` package) for the given conditions and
return the neutral atmosphere as a [`NeutralProfile`](@ref). Unlike [`find_msis_file`](@ref),
nothing is written to disk — the computed profile lives in the returned struct and round-trips
through `physics_state.jld2`.

This is the neutral counterpart of [`run_iri`](@ref). Levels where the model does not report a
species (it returns `NaN` for N and anomalous O at low altitude) are dropped for that species
only, so each density keeps just the altitudes where it is defined.

# Example
```julia
neutrals = run_msis(; year=2005, month=10, day=8, hour=22, minute=0, lat=69.58, lon=19.23)
electrons = run_iri(; year=2005, month=10, day=8, hour=22, minute=0, lat=69.58, lon=19.23)
model    = AuroraModel(altitude_lims, θ_lims, E_max, neutrals, electrons)
```
"""
function run_msis(; year = 2018, month = 12, day = 7, hour = 11, minute = 15,
                  lat = 76, lon = 5, height = 85:1:700, verbose = true)
    msis_data, _ = calculate_msis_data(; year, month, day, hour, minute, lat, lon, height,
                                        verbose)
    data    = msis_data[2:end, :]                  # drop the header row
    h_m     = Float64.(data[:, 1]) .* 1e3          # height(km) → m
    instant = DateTime(year, month, day, hour, minute)
    label   = "NRLMSIS 2.1 $instant $(lat)N/$(lon)E"

    # Columns of calculate_msis_data, all number densities in m⁻³
    densities = (:N2 => 3, :O2 => 4, :O => 5, :He => 6, :H => 7, :Ar => 8, :N => 9,
                 :NO => 11)
    columns   = Dict(s => Float64.(data[:, c]) for (s, c) in densities)

    densities, dropped = species_densities(h_m, columns, label)
    return NeutralProfile(densities; source = label, dropped)
end

"""
    read_msis_file(msis_file) -> NeutralProfile

Read every species from an MSIS text file generated by AURORA and return them as a
[`NeutralProfile`](@ref) on the file's native altitude grid. The file is read once, here, so
the result is self-contained and round-trips with no path dependency.

Kept for backward compatibility with existing MSIS files; for new runs prefer
[`run_msis`](@ref). Index the result to reach one species: `read_msis_file(f)[:N2]`.

# Example
```julia
neutrals = read_msis_file(msis_file)
model    = AuroraModel(altitude_lims, θ_lims, E_max, neutrals, electrons)
```
"""
function read_msis_file(msis_file::AbstractString)
    raw   = load_msis(msis_file)
    h_m   = raw.data.height_km .* 1e3
    label = "MSIS file $(basename(msis_file))"

    columns = Dict(species => getproperty(raw.data, species)
                   for species in (:N2, :O2, :O, :He, :H, :Ar, :N, :NO)
                   if hasproperty(raw.data, species))
    densities, dropped = species_densities(h_m, columns, label)
    return NeutralProfile(densities; source = label, dropped)
end

"""
    read_ccmc_msis(file) -> NeutralProfile

Read the neutral atmosphere from a CCMC ModelWeb NRLMSIS text export and return it as a
[`NeutralProfile`](@ref). The CCMC table has a variable-length preamble, a single-line column
header, densities in cm⁻³, and a `9.999E-38` sentinel for species that the model does not
report at a given altitude; this reader locates the header (the line containing `N2den`),
converts cm⁻³ → m⁻³, and drops sentinel levels per species.

Columns are resolved by their header name (`Heit(km)`, `N2den(cm-3)`, …) rather than by
position, so any export using these names is read correctly whatever its column order, and a
species absent from the export is simply not returned. This covers the NRLMSIS 2.x and
NRLMSISE-00 exports, which use the same column names.

# Example
```julia
neutrals = read_ccmc_msis("nrlmsis_output.txt")
model    = AuroraModel(altitude_lims, θ_lims, E_max, neutrals, electrons)
```
"""
function read_ccmc_msis(file::AbstractString)
    lines      = readlines(file)
    header_idx = findfirst(l -> occursin("N2den", l), lines)
    header_idx === nothing && throw(ArgumentError(
        "read_ccmc_msis: could not find the CCMC column header (a line containing " *
        "\"N2den\") in $file. Is this a CCMC ModelWeb NRLMSIS export?"))

    # The header names its columns one-for-one with the data columns, so look up the ones we
    # need by name. The unit is part of the name, which keeps a change of unit from being
    # read as if it were cm⁻³.
    header = split(lines[header_idx])
    column = Dict(name => i for (i, name) in enumerate(header))
    columns_found = "Columns found: " * join(header, ", ")

    haskey(column, "Heit(km)") || throw(ArgumentError(
        "read_ccmc_msis: no altitude column \"Heit(km)\" in the header of $file.\n" *
        columns_found))
    h_col = column["Heit(km)"]

    species_columns = [s => column[name] for (s, name) in
                       (:O  => "Oden(cm-3)",  :N2 => "N2den(cm-3)", :O2 => "O2den(cm-3)",
                        :NO => "NOden(cm-3)", :He => "Heden(cm-3)", :Ar => "Arden(cm-3)",
                        :H  => "Hden(cm-3)",  :N  => "Nden(cm-3)") if haskey(column, name)]
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
        # An unparseable field costs that species this level, not the level for every
        # species — the same per-species rule the sentinel handling below follows.
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
    return NeutralProfile(densities; source = label, dropped)
end
