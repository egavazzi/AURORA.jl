# ======================================================================================== #
#                              CollisionChannel                                            #
# ======================================================================================== #

"""
    CollisionChannel{F}

One inelastic collision channel of a neutral species: a cross section, the energy it removes
from the colliding electron, and the number of secondary electrons it ejects.

In a species, `channels[i]` is row `i + 1` of `cross_sections` and `excitation_levels` (row
1 is elastic).

# Fields
- `name::String`: channel identifier, e.g. `"a3sup"`, `"ionx2sgp"`. Used for labels and for
    [`channel`](@ref) lookup.
- `cross_section::F`: callable mapping a vector of energies (eV) to cross sections (m²).
    Must be a named function, a functor, or an [`ExprLaw`](@ref).
- `energy_loss::Float64`: energy removed from the electron per collision (eV).
- `n_secondaries::Int`: secondary electrons ejected — `0` for excitation, `1` for single
    ionization, `2` for double ionization.
- `source::String`: free-form provenance (reference, table, digitization note). May be empty.

    CollisionChannel(name, cross_section, energy_loss, n_secondaries; source = "")

Build a channel. `energy_loss` must be non-negative, and strictly positive for an ionizing
channel.

    CollisionChannel(c::CollisionChannel; name, cross_section, energy_loss, n_secondaries, source)

Copy `c`, replacing the given fields.

# Example
```julia
CollisionChannel("ionx2sgp", e_N2ionx2sgp, 15.581, 1)
CollisionChannel("mystate", @law(E -> 1e-21 .* (E .> 5)), 5.0, 0; source = "made up")
```
"""
struct CollisionChannel{F}
    name::String
    cross_section::F
    energy_loss::Float64
    n_secondaries::Int
    source::String

    function CollisionChannel{F}(name, cross_section, energy_loss, n_secondaries,
                                 source) where {F}
        channel_name = String(name)
        require_reproducible(cross_section, "cross_section of channel $(channel_name)")
        loss = Float64(energy_loss)
        isfinite(loss) && loss >= 0 || throw(ArgumentError(
            "channel $(channel_name) has an energy loss of $(loss) eV; it must be finite and \
             non-negative"))
        n_sec = Int(n_secondaries)
        0 <= n_sec <= 2 || throw(ArgumentError(
            "channel $(channel_name) ejects $(n_sec) secondary electrons; a collision \
             channel must eject 0 (excitation), 1 (single ionization) or 2 (double \
             ionization)"))
        n_sec == 0 || loss > 0 || throw(ArgumentError(
            "channel $(channel_name) ejects $(n_sec) secondary electrons but costs no \
             energy; an ionizing channel needs a positive energy loss"))
        return new{F}(channel_name, cross_section, loss, n_sec, String(source))
    end
end

function CollisionChannel(name, cross_section, energy_loss, n_secondaries; source = "")
    return CollisionChannel{typeof(cross_section)}(name, cross_section, energy_loss,
                                                   n_secondaries, source)
end

function CollisionChannel(c::CollisionChannel;
                          name = c.name,
                          cross_section = c.cross_section,
                          energy_loss = c.energy_loss,
                          n_secondaries = c.n_secondaries,
                          source = c.source)
    return CollisionChannel(name, cross_section, energy_loss, n_secondaries; source)
end

function Base.show(io::IO, c::CollisionChannel)
    print(io, "CollisionChannel(\"", c.name, "\", ", c.energy_loss, " eV, n_secondaries=",
          c.n_secondaries, ")")
end

"""
    channel_names(channels) → Vector{String}
    channel_names(sp::NeutralSpecies) → Vector{String}

Names of the inelastic collision channels, in row order.
"""
channel_names(channels) = [c.name for c in channels]

# Channel names are lookup keys and output labels, so a species must not repeat one.
function check_channel_names(channels)
    names = channel_names(channels)
    allunique(names) && return nothing
    duplicated = unique(filter(n -> count(==(n), names) > 1, names))
    throw(ArgumentError("channel names must be distinct; duplicated: $(duplicated)"))
end

"""
    ionizing_channels(channels)
    ionizing_channels(sp::NeutralSpecies)

The channels that eject at least one secondary electron, in row order.
"""
function ionizing_channels(channels)
    return filter(c -> c.n_secondaries > 0, channels)
end

"""
    channel(channels, name::AbstractString) → CollisionChannel
    channel(sp::NeutralSpecies, name::AbstractString) → CollisionChannel

The channel called `name`. Throws a `KeyError` when there is none, and an `ArgumentError`
when several channels share the name.
"""
function channel(channels, name::AbstractString)
    i = findall(c -> c.name == name, channels)
    isempty(i) && throw(KeyError(name))
    length(i) == 1 || throw(ArgumentError("Multiple channels are named $(name)"))
    return channels[only(i)]
end

"""
    channel_cross_sections(elastic_cross_section, channels, E_centers) → Matrix{Float64}

Build the `[n_levels × n_E]` cross-section matrix of a species: row 1 is the elastic cross
section, row `i + 1` is `channels[i]`.
"""
function channel_cross_sections(elastic_cross_section, channels, E_centers::AbstractVector)
    Base.require_one_based_indexing(E_centers, channels)
    σ = zeros(length(channels) + 1, length(E_centers))
    σ[1, :] .= elastic_cross_section(E_centers)
    for (i, c) in enumerate(channels)
        σ[i + 1, :] .= c.cross_section(E_centers)
    end
    return σ
end

"""
    channel_excitation_levels(channels) → Matrix{Float64}

Build the `[n_levels × 2]` excitation-level matrix of a species: column 1 is the energy loss
(eV), column 2 the number of secondary electrons. Row 1 is the elastic channel (zero, zero),
row `i + 1` is `channels[i]`.
"""
function channel_excitation_levels(channels)
    Base.require_one_based_indexing(channels)
    levels = zeros(length(channels) + 1, 2)
    for (i, c) in enumerate(channels)
        levels[i + 1, 1] = c.energy_loss
        levels[i + 1, 2] = c.n_secondaries
    end
    return levels
end
