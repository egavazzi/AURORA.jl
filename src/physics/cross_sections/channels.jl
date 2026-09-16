include("channels_N2.jl")
include("channels_O2.jl")
include("channels_O.jl")

"""
    default_channels(species::Symbol) → Vector{CollisionChannel}

The built-in inelastic collision channels of `:N2`, `:O2` or `:O`, in the row order of the
`cross_sections` and `excitation_levels` matrices.

The tables live in `src/physics/cross_sections/channels_N2.jl`, `channels_O2.jl` and
`channels_O.jl`.
"""
function default_channels(species::Symbol)
    species === :N2 && return default_channels_N2()
    species === :O2 && return default_channels_O2()
    species === :O  && return default_channels_O()
    throw(ArgumentError("No built-in collision channels for $(species); AURORA ships tables \
                         for :N2, :O2 and :O. Build a `Vector{CollisionChannel}` to describe \
                         another species."))
end

"""
    default_elastic_cross_section(species::Symbol)

The built-in elastic cross section of `:N2`, `:O2` or `:O`.
"""
function default_elastic_cross_section(species::Symbol)
    species === :N2 && return default_elastic_cross_section_N2()
    species === :O2 && return default_elastic_cross_section_O2()
    species === :O  && return default_elastic_cross_section_O()
    throw(ArgumentError("No built-in elastic cross section for $(species); AURORA ships \
                         data for :N2, :O2 and :O"))
end

"""
    default_secondary_law(species::Symbol)

The built-in secondary-electron energy distribution of `:N2`, `:O2` or `:O`, for building the
cascading transfer matrices.
"""
function default_secondary_law(species::Symbol)
    species === :N2 && return default_secondary_law_N2()
    species === :O2 && return default_secondary_law_O2()
    species === :O  && return default_secondary_law_O()
    throw(ArgumentError("No built-in secondary-electron law for $(species); AURORA ships \
                         laws for :N2, :O2 and :O"))
end

"""
    default_cascading_spec(species::Symbol) → CascadingSpec

The cascading spec of `:N2`, `:O2` or `:O`, derived from the built-in channel table and
secondary-electron law.
"""
function default_cascading_spec(species::Symbol)
    return CascadingSpec(String(species), default_secondary_law(species);
                         channels = default_channels(species))
end
