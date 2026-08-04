"""
    Ionosphere{FT, V<:AbstractVector{FT}, S}

Electron background: electron density and temperature, sampled on the model altitude grid from
an `electron_source` (an [`ElectronProfile`](@ref) or any callable `h_atm (m) → (ne, Te)`).
Neutral species densities are owned by the individual [`NeutralSpecies`](@ref) objects.
"""
struct Ionosphere{FT, V<:AbstractVector{FT}, S}
    Te::V
    ne::V
    atmosphere_source::String
    electron_source::S
end

function Ionosphere(atmosphere_source::AbstractString, electron_source,
                    h_atm::AbstractVector)
    es     = to_electron_source(electron_source)   # legacy file path → ElectronProfile
    ne, Te = es(h_atm)
    FT     = eltype(Te)
    return Ionosphere{FT, typeof(Te), typeof(es)}(Te, ne, string(atmosphere_source), es)
end

function Base.show(io::IO, iono::Ionosphere)
    print(io, "Ionosphere($(length(iono.ne)) altitudes)")
end

function Base.show(io::IO, ::MIME"text/plain", iono::Ionosphere)
    nz = length(iono.ne)
    println(io, "Ionosphere:")
    println(io, "├── Altitudes: $(nz)")
    println(io, "├── Atmosphere: $(iono.atmosphere_source)")
    println(io, "├── Electrons: $(electron_source_label(iono.electron_source))")
    println(io, "├── Max Te:    $(round(maximum(iono.Te), sigdigits=3)) K")
    print(io,   "└── Max ne:    $(round(maximum(iono.ne), sigdigits=3)) m⁻³")
end
