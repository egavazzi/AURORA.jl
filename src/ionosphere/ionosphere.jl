"""
    Ionosphere{FT, V<:AbstractVector{FT}}

Electron and thermal background: electron density and temperature plus neutral temperature.
Neutral species densities are owned by the individual [`NeutralSpecies`](@ref) objects.
"""
struct Ionosphere{FT, V<:AbstractVector{FT}}
    Tn::V
    Te::V
    ne::V
    msis_file::String
    iri_file::String
end

function Ionosphere(msis_file::AbstractString, iri_file::AbstractString,
                    h_atm::AbstractVector)
    msis_raw = load_msis(msis_file)
    msis = interpolate_msis_to_grid(msis_raw.data, h_atm)
    Tn = msis.T
    ne, Te = load_electron_densities(iri_file, h_atm)
    FT = eltype(Tn)
    return Ionosphere{FT, typeof(Tn)}(Tn, Te, ne, string(msis_file), string(iri_file))
end

"""
    load_electron_densities(iri_file, h_atm)

Load the electron density and temperature from an IRI file that was generated and saved
using AURORA's IRI interface. Then interpolate the profiles over AURORA's altitude grid.

# Calling
`ne, Te = load_electron_densities(iri_file, h_atm)`

# Inputs
- `iri_file`: absolute path to the iri file to read ne and Te from. String
- `h_atm`: altitude (m). Vector [nZ]

# Returns
- `ne`: e- density (m⁻³). Vector [nZ]
- `Te`: e- temperature (K). Vector [nZ]

# See also
[`load_iri`](@ref), [`interpolate_iri_to_grid`](@ref)
"""
function load_electron_densities(iri_file, h_atm)
    iri_raw = load_iri(iri_file)
    iri = interpolate_iri_to_grid(iri_raw.data, h_atm)

    return iri.ne, iri.Te
end

function Base.show(io::IO, iono::Ionosphere)
    print(io, "Ionosphere($(length(iono.Tn)) altitudes)")
end

function Base.show(io::IO, ::MIME"text/plain", iono::Ionosphere)
    nz = length(iono.Tn)
    println(io, "Ionosphere:")
    println(io, "├── Altitudes: $(nz)")
    println(io, "├── MSIS file: $(basename(iono.msis_file))")
    println(io, "├── IRI  file: $(basename(iono.iri_file))")
    println(io, "├── Max Tn:    $(round(maximum(iono.Tn), sigdigits=3)) K")
    println(io, "├── Max Te:    $(round(maximum(iono.Te), sigdigits=3)) K")
    print(io,   "└── Max ne:    $(round(maximum(iono.ne), sigdigits=3)) m⁻³")
end
