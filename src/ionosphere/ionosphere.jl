using SpecialFunctions: erf

"""
    Ionosphere{FT, V<:AbstractVector{FT}}

Atmospheric state containing neutral and electron densities and temperatures.
"""
struct Ionosphere{FT, V<:AbstractVector{FT}}
    nN2::V
    nO2::V
    nO::V
    Tn::V
    Te::V
    ne::V
    msis_file::String
    iri_file::String
end

function Ionosphere(msis_file::AbstractString, iri_file::AbstractString,
                    h_atm::AbstractVector)
    n_neutrals, Tn = load_neutral_densities(msis_file, h_atm)
    ne, Te = load_electron_densities(iri_file, h_atm)
    FT = eltype(Tn)
    return Ionosphere{FT, typeof(Tn)}(
        n_neutrals.nN2, n_neutrals.nO2, n_neutrals.nO,
        Tn, Te, ne,
        string(msis_file), string(iri_file)
    )
end

"""
    n_neutrals(iono::Ionosphere)

Return neutral densities as a NamedTuple `(; nN2, nO2, nO)`.
"""
n_neutrals(iono::Ionosphere) = (; nN2=iono.nN2, nO2=iono.nO2, nO=iono.nO)

"""
    load_neutral_densities(msis_file, h_atm)

Load the neutral densities and temperature from a MSIS file that was generated and saved
using AURORA's MSIS interface. Then interpolate the profiles over AURORA's altitude grid.

Upper boundary conditions are applied to smoothly transition the densities to zero.

# Calling
`n_neutrals, Tn = load_neutral_densities(msis_file, h_atm)`

# Inputs
- `msis_file`: absolute path to the msis file to read n_neutrals and Tn from. String
- `h_atm`: altitude (m). Vector [nZ]

# Returns
- `n_neutrals`: neutral densities (m⁻³). Named tuple of vectors ([nZ], ..., [nZ])
- `Tn`: neutral temperature (K). Vector [nZ]

# See also
[`load_msis`](@ref), [`interpolate_msis_to_grid`](@ref)
"""
function load_neutral_densities(msis_file, h_atm)
    msis_raw = load_msis(msis_file)
    msis = interpolate_msis_to_grid(msis_raw.data, h_atm)
    nN2 = msis.N2
    nO2 = msis.O2
    nO = msis.O

    # Apply upper boundary conditions to smoothly transition densities to zero
    # Set the last 3 points to zero
    nN2[end-2:end] .= 0
    nO2[end-2:end] .= 0
    nO[end-2:end] .= 0

    # Apply smooth transition using error function for the 3 points below
    erf_factor = (erf.((1:-1:-1) / 2) .+ 1) / 2
    nN2[end-5:end-3] .= erf_factor .* nN2[end-5:end-3]
    nO2[end-5:end-3] .= erf_factor .* nO2[end-5:end-3]
    nO[end-5:end-3] .= erf_factor .* nO[end-5:end-3]

    # Ensure no negative densities
    if any(nN2 .< 0) || any(nO2 .< 0) || any(nO .< 0)
        @warn "Negative densities found. They were set to 0, but you might want to check " *
              "why that happened"
        nN2[nN2 .< 0] .= 0
        nO2[nO2 .< 0] .= 0
        nO[nO .< 0] .= 0
    end

    n_neutrals = (; nN2, nO2, nO)
    return n_neutrals, msis.T
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
    print(io, "Ionosphere($(length(iono.nN2)) altitudes)")
end

function Base.show(io::IO, ::MIME"text/plain", iono::Ionosphere)
    nz = length(iono.nN2)
    println(io, "Ionosphere:")
    println(io, "├── Altitudes: $(nz)")
    println(io, "├── MSIS file: $(basename(iono.msis_file))")
    println(io, "├── IRI  file: $(basename(iono.iri_file))")
    println(io, "├── Max n(N2): $(round(maximum(iono.nN2), sigdigits=3)) m^-3")
    println(io, "├── Max n(O2): $(round(maximum(iono.nO2), sigdigits=3)) m^-3")
    print(io, "└── Max n(O):  $(round(maximum(iono.nO), sigdigits=3)) m^-3")
end
