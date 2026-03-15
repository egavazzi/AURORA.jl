"""
    CrossSectionData{NT1<:NamedTuple, NT2<:NamedTuple}

Collision cross-sections and excitation thresholds for all neutral species.
"""
struct CrossSectionData{NT1<:NamedTuple, NT2<:NamedTuple}  # More explicit
    σ_neutrals::NT1
    collision_levels::NT2
end

function CrossSectionData(energy_grid::EnergyGrid)
    σ_neutrals = load_cross_sections(energy_grid)
    collision_levels = load_excitation_threshold()
    return CrossSectionData(σ_neutrals, collision_levels)
end

function Base.show(io::IO, cs::CrossSectionData)
    nE = size(cs.σ_neutrals.σ_N2, 2)
    print(io, "CrossSectionData($(nE) energies)")
end

function Base.show(io::IO, ::MIME"text/plain", cs::CrossSectionData)
    nE = size(cs.σ_neutrals.σ_N2, 2)
    n_levels_N2 = size(cs.collision_levels.N2_levels, 1)
    n_levels_O2 = size(cs.collision_levels.O2_levels, 1)
    n_levels_O = size(cs.collision_levels.O_levels, 1)
    println(io, "CrossSectionData:")
    println(io, "├── Energy bins: $(nE)")
    println(io, "├── N2 levels: $(n_levels_N2)")
    println(io, "├── O2 levels: $(n_levels_O2)")
    print(io, "└── O levels: $(n_levels_O)")
end
