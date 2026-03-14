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
