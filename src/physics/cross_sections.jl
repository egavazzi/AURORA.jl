"""
    CrossSectionData{NT1<:NamedTuple, NT2<:NamedTuple}

Collision cross-sections and excitation thresholds for all neutral species.
"""
struct CrossSectionData{NT1<:NamedTuple, NT2<:NamedTuple}  # More explicit
    σ_neutrals::NT1
    E_levels_neutrals::NT2
end

function CrossSectionData(E::AbstractVector, dE::AbstractVector)
    σ_neutrals = load_cross_sections(E, dE)
    E_levels_neutrals = load_excitation_threshold()
    return CrossSectionData(σ_neutrals, E_levels_neutrals)
end

function Base.show(io::IO, cs::CrossSectionData)
    nE = size(cs.σ_neutrals.σ_N2, 2)
    print(io, "CrossSectionData($(nE) energies)")
end

function Base.show(io::IO, ::MIME"text/plain", cs::CrossSectionData)
    nE = size(cs.σ_neutrals.σ_N2, 2)
    n_levels_N2 = size(cs.E_levels_neutrals.N2_levels, 1)
    n_levels_O2 = size(cs.E_levels_neutrals.O2_levels, 1)
    n_levels_O = size(cs.E_levels_neutrals.O_levels, 1)
    println(io, "CrossSectionData:")
    println(io, "├── Energy bins: $(nE)")
    println(io, "├── N2 levels: $(n_levels_N2)")
    println(io, "├── O2 levels: $(n_levels_O2)")
    print(io, "└── O levels: $(n_levels_O)")
end
