include("channels_N2.jl")
include("channels_O2.jl")
include("channels_O.jl")

"""
    default_collisions(species::Symbol) → NamedTuple

The built-in collision data of `:N2`, `:O2` or `:O`, as the keyword arguments
`(; elastic_cross_section, channels, secondary_law)` of [`NeutralSpecies`](@ref). The
channels are in the row order of the `cross_sections` and `excitation_levels` matrices.
"""
function default_collisions(species::Symbol)
    species === :N2 && return default_collisions_N2()
    species === :O2 && return default_collisions_O2()
    species === :O  && return default_collisions_O()
    throw(ArgumentError("No built-in collision channels for $(species); AURORA ships tables \
                         for :N2, :O2 and :O. Build a `Vector{CollisionChannel}` to describe \
                         another species."))
end

"""
    default_cascading_spec(species::Symbol) → CascadingSpec

The cascading spec of `:N2`, `:O2` or `:O`, derived from the built-in channel table and
secondary-electron law.
"""
function default_cascading_spec(species::Symbol)
    d = default_collisions(species)
    return cascading_spec_from_channels(String(species), d.secondary_law, d.channels)
end
