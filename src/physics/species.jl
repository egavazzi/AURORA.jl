# ======================================================================================== #
#                              NeutralSpecies                                              #
# ======================================================================================== #

"""
    NeutralSpecies

All per-species data needed to advance the transport equation through one neutral species.

The collision physics is described by `elastic_cross_section`, `channels` and
`secondary_law`. Everything the solvers read — `cross_sections`, `excitation_levels`,
`cascading_spec` and `cascading_data` — is derived from those three at `initialize!(model)`,
so a channel added, removed or edited before that call is reflected everywhere at once.

# Fields
- `name::Symbol`: short identifier (e.g. `:N2`, `:O2`, `:O`)
- `density_source`: callable `h_atm (m) → density (m⁻³)` that (re)samples `density`.
    Can be a [`DensityProfile`](@ref) or any callable.
    Untyped so it can be replaced freely before calling `initialize!(model)`.
- `density::Vector{Float64}`: density profile sampled on the model altitude grid (m⁻³).
    Empty until `initialize!(model)` is called.
- `elastic_cross_section`: callable mapping energies (eV) to the elastic cross section (m²).
- `channels::Vector{CollisionChannel}`: the inelastic collision channels, in the row order
    of `cross_sections` and `excitation_levels`.
- `secondary_law`: callable `(E_secondary, E_primary) -> Float64` giving the
    secondary-electron energy distribution of the ionizing channels.
- `phase_fcn_generator`: callable `(θ, E) -> (phaseE, phaseI)` that (re)builds
    `phase_fcn` whenever the pitch-angle or energy grid changes.
    Untyped so it can be replaced freely before calling `initialize!(model)`.
- `phase_fcn`: tuple `(phaseE, phaseI)` of `[n_θ × n_E]` matrices materialized from the
    generator on the model's scattering θ grid and energy centers.
    Holds 0×0 placeholder matrices until `initialize!(model)` is called.
- `cross_sections::Matrix{Float64}`: cross sections evaluated on the model energy grid,
    shape `[n_levels × n_E]` (m²). Row 1 is elastic, row `i + 1` is `channels[i]`.
    Empty until `initialize!(model)` is called.
- `excitation_levels::Matrix{Float64}`: energy losses and secondary counts, shape
    `[n_levels × 2]`, in the same row order as `cross_sections`.
    Empty until `initialize!(model)` is called.
- `cascading_spec::CascadingSpec`: ionization thresholds and secondary-distribution law the
    cascading transfer matrices are built from.
- `cascading_data::SpeciesCascadingCache`: cascading transfer matrices, filled in by
    `load_or_compute_cascading!`

!!! note "Mutation happens before `initialize!`"
    Assigning `channels`, `elastic_cross_section`, `secondary_law` or `density_source` does
    not mark the model uninitialized, because a species does not know which model holds it.
    Edit a species of an already-initialized model, then call `initialize!(model)` explicitly
    to rebuild the derived data.
"""
mutable struct NeutralSpecies
    name::Symbol
    density_source         # untyped: freely replaceable before initialize!(model)
    density::Vector{Float64}
    elastic_cross_section  # untyped: freely replaceable before initialize!(model)
    channels::Vector{CollisionChannel}
    secondary_law          # untyped: freely replaceable before initialize!(model)
    phase_fcn_generator    # untyped: freely replaceable before initialize!(model)
    phase_fcn::Tuple{Matrix{Float64}, Matrix{Float64}}
    cross_sections::Matrix{Float64}
    excitation_levels::Matrix{Float64}
    cascading_spec::CascadingSpec
    cascading_data::SpeciesCascadingCache

    function NeutralSpecies(name, density_source, density, elastic_cross_section, channels,
                            secondary_law, phase_fcn_generator, phase_fcn, cross_sections,
                            excitation_levels, cascading_spec, cascading_data)
        require_reproducible(density_source, "density_source")
        require_reproducible(elastic_cross_section, "elastic_cross_section")
        require_reproducible(secondary_law, "secondary_law")
        require_reproducible(phase_fcn_generator, "phase_fcn_generator")
        return new(Symbol(name), density_source, density, elastic_cross_section,
                   collect(CollisionChannel, channels), secondary_law,
                   phase_fcn_generator, phase_fcn, cross_sections, excitation_levels,
                   cascading_spec, cascading_data)
    end
end

"""
    NeutralSpecies(name::Symbol, density_source; elastic_cross_section, channels,
                   secondary_law, phase_fcn_generator)

Build a lightweight `NeutralSpecies`. The grid-dependent fields (`density`,
`cross_sections`, `excitation_levels`, `phase_fcn`) start empty and are populated by
`initialize!(model)`.

# Example
```julia
gas = NeutralSpecies(:MyGas, @law(z -> 1e18 .* exp.(-z ./ 30e3));
                     elastic_cross_section = AURORA.e_N2elastic,
                     channels = [CollisionChannel("exc", AURORA.e_N2a3sup, 6.17, 0),
                                 CollisionChannel("ion", AURORA.e_N2ionx2sgp, 15.6, 1)],
                     secondary_law = @law((E_s, E_p) -> 1.0 / (11.4^2 + E_s^2)),
                     phase_fcn_generator = phase_fcn_N2)
```
"""
function NeutralSpecies(name::Symbol, density_source; elastic_cross_section, channels,
                        secondary_law, phase_fcn_generator)
    empty_mat = Matrix{Float64}(undef, 0, 0)
    spec = CascadingSpec(String(name), secondary_law; channels)
    return NeutralSpecies(
        name,
        density_source,
        Float64[],
        elastic_cross_section,
        channels,
        secondary_law,
        phase_fcn_generator,
        (empty_mat, copy(empty_mat)),
        copy(empty_mat),
        copy(empty_mat),
        spec,
        SpeciesCascadingCache(spec),
    )
end

# Density profiles, cross sections and phase-function generators are commonly swapped in via
# direct field assignment, which bypasses the constructor. Intercept those assignments to
# enforce the reproducibility rule there too.
function Base.setproperty!(sp::NeutralSpecies, name::Symbol, value)
    if name in (:density_source, :elastic_cross_section, :phase_fcn_generator)
        require_reproducible(value, String(name))
    elseif name === :secondary_law
        require_reproducible(value, "secondary_law")
    end
    ty = fieldtype(typeof(sp), name)
    return setfield!(sp, name, value isa ty ? value : convert(ty, value))
end

function Base.getindex(species::Tuple{Vararg{NeutralSpecies}}, name::Symbol)
    found_index = 0
    for (i, sp) in pairs(species)
        if sp.name == name
            found_index == 0 || throw(ArgumentError("Multiple species are named $(name)"))
            found_index = i
        end
    end
    found_index == 0 && throw(KeyError(name))
    return species[found_index]
end

channel_names(sp::NeutralSpecies) = channel_names(sp.channels)
ionizing_channels(sp::NeutralSpecies) = ionizing_channels(sp.channels)
channel(sp::NeutralSpecies, name::AbstractString) = channel(sp.channels, name)

"""
    rebuild_collision_data!(sp::NeutralSpecies, E_centers)

Derive `cross_sections`, `excitation_levels` and the cascading spec of `sp` from its channel
table, evaluating the cross sections on `E_centers`.

A spec equal to the one the cascading cache already holds leaves that cache object in place,
for `load_or_compute_cascading!` to refresh against the current energy grid; a spec that
differs replaces the cache with an empty one built from it.
"""
function rebuild_collision_data!(sp::NeutralSpecies, E_centers::AbstractVector)
    sp.cross_sections    = channel_cross_sections(sp.elastic_cross_section, sp.channels,
                                                  E_centers)
    sp.excitation_levels = channel_excitation_levels(sp.channels)
    spec = CascadingSpec(String(sp.name), sp.secondary_law; sp.channels)
    if !describes_same_cascading(sp.cascading_spec, spec)
        sp.cascading_spec = spec
        sp.cascading_data = SpeciesCascadingCache(spec)
    end
    return nothing
end

# Two specs build the same transfer matrices when they agree on the thresholds, the secondary
# counts and the law object itself.
function describes_same_cascading(a::CascadingSpec, b::CascadingSpec)
    return a.name == b.name &&
           a.ionization_thresholds == b.ionization_thresholds &&
           a.n_secondaries == b.n_secondaries &&
           a.secondary_law === b.secondary_law
end


# ======================================================================================== #
#                        Default species convenience constructors                          #
# ======================================================================================== #

"""
    N2Species(density_source)
    N2Species(neutrals::NeutralAtmosphere)
    N2Species(msis_file::AbstractString)

Construct the default N₂ species: elastic cross section, collision channels and
secondary-electron law from [`default_channels`](@ref) and friends, phase function from
[`phase_fcn_N2`](@ref).

`density_source` can be a [`DensityProfile`](@ref) or any callable `h_atm (m) → density
(m⁻³)`. Passing a [`NeutralAtmosphere`](@ref) is shorthand for `neutrals[:N2]`. Passing an MSIS
file path string is shorthand for `read_msis_file(msis_file)[:N2]`. Grid-dependent fields are
populated later by `initialize!(model)`.
"""
function N2Species(neutrals::NeutralAtmosphere)
    return N2Species(neutrals[:N2])
end
function N2Species(msis_file::AbstractString)
    return N2Species(read_msis_file(msis_file)[:N2])
end
function N2Species(density_source)
    return default_species(:N2, density_source, phase_fcn_N2)
end

"""
    O2Species(density_source)
    O2Species(neutrals::NeutralAtmosphere)
    O2Species(msis_file::AbstractString)

Default O₂ species, analogous to [`N2Species`](@ref).
"""
function O2Species(neutrals::NeutralAtmosphere)
    return O2Species(neutrals[:O2])
end
function O2Species(msis_file::AbstractString)
    return O2Species(read_msis_file(msis_file)[:O2])
end
function O2Species(density_source)
    return default_species(:O2, density_source, phase_fcn_O2)
end

"""
    OSpecies(density_source)
    OSpecies(neutrals::NeutralAtmosphere)
    OSpecies(msis_file::AbstractString)

Default O species, analogous to [`N2Species`](@ref).
"""
function OSpecies(neutrals::NeutralAtmosphere)
    return OSpecies(neutrals[:O])
end
function OSpecies(msis_file::AbstractString)
    return OSpecies(read_msis_file(msis_file)[:O])
end
function OSpecies(density_source)
    return default_species(:O, density_source, phase_fcn_O)
end

function default_species(name::Symbol, density_source, phase_fcn_generator)
    return NeutralSpecies(name, density_source;
                          elastic_cross_section = default_elastic_cross_section(name),
                          channels              = default_channels(name),
                          secondary_law         = default_secondary_law(name),
                          phase_fcn_generator)
end


# ======================================================================================== #
#                              Display                                                     #
# ======================================================================================== #

function Base.show(io::IO, sp::NeutralSpecies)
    print(io, "NeutralSpecies($(sp.name), ", length(sp.density), " altitudes)")
end

function Base.show(io::IO, ::MIME"text/plain", sp::NeutralSpecies)
    println(io, "NeutralSpecies(", sp.name, "):")
    println(io, "├── Density source:  ", profile_label(sp.density_source))
    println(io, "├── Channels:        ", length(sp.channels), " (",
                                         length(ionizing_channels(sp)), " ionizing)")
    if isempty(sp.density)
        println(io, "├── (not initialized — call initialize!(model))")
    else
        println(io, "├── Altitudes:        ", length(sp.density))
        println(io, "├── Max density:      ", round(maximum(sp.density), sigdigits=3), " m⁻³")
        println(io, "├── Cross sections:   ", size(sp.cross_sections, 1), " levels × ",
                                              size(sp.cross_sections, 2), " energies")
    end
    print(io,   "└── Cascading:        ", sp.cascading_spec.name,
                " (", length(sp.cascading_spec.ionization_thresholds), " thresholds)")
end
