using SpecialFunctions: erf
using DataInterpolations: PCHIPInterpolation, ExtrapolationType

# ======================================================================================== #
#                              Density-profile types                                       #
# ======================================================================================== #

"""
    MSISDensity(msis_file, species)

Density profile backed by an MSIS file. Callable on an altitude grid (m), returns the
density vector (m⁻³) for the requested species (`:N2`, `:O2`, or `:O`).

# Example
```julia
profile = MSISDensity(find_msis_file(), :N2)
n = profile(altitude_grid.h)   # Vector{Float64}, length == length(altitude_grid.h)
```
"""
struct MSISDensity
    msis_file::String
    species::Symbol
end

(d::MSISDensity)(h_atm::AbstractVector) = load_msis_density(d.msis_file, d.species, h_atm)

"""
    VectorDensity(h, n)

Density profile defined by user-supplied altitude (`h`, m) and density (`n`, m⁻³) vectors.
Callable on any altitude grid (m); evaluates via PCHIP interpolation in log-space,
consistent with AURORA's MSIS interpolation convention.

# Example
```julia
profile = VectorDensity(h_msis_m, n_N2)
n = profile(altitude_grid.h)
```
"""
struct VectorDensity
    h::Vector{Float64}   # altitude (m)
    n::Vector{Float64}   # density (m⁻³)
end

function (d::VectorDensity)(h_atm::AbstractVector)
    itp = PCHIPInterpolation(log.(d.n), d.h ./ 1e3;
                             extrapolation = ExtrapolationType.Linear)
    return exp.(itp(collect(Float64, h_atm) ./ 1e3))
end


# ======================================================================================== #
#                              NeutralSpecies                                              #
# ======================================================================================== #

"""
    NeutralSpecies{F, S<:CascadingSpec, C<:SpeciesCascadingCache, P}

All per-species data needed to advance the transport equation through one neutral species.

# Fields
- `name::Symbol`: short identifier (e.g. `:N2`, `:O2`, `:O`)
- `density_profile`: callable `h_atm (m) → density (m⁻³)` used to (re)sample `density`.
    Can be an [`MSISDensity`](@ref), a [`VectorDensity`](@ref), or any callable.
    This field is untyped so it can be replaced freely before calling `initialize!(model)`.
- `density::Vector{Float64}`: density profile sampled on the model altitude grid (m⁻³),
    with an erf-tail boundary applied at the top. Empty until `initialize!(model)` is called.
- `cross_sections::Matrix{Float64}`: collision cross sections, shape `[n_levels × n_E]` (m²).
    Empty until `initialize!(model)` is called.
- `excitation_levels::Matrix{Float64}`: excitation thresholds and secondary counts,
    shape `[n_levels × 2]`. First column = threshold energy (eV), second column = number
    of secondaries produced (non-zero only for ionizing channels).
    Empty until `initialize!(model)` is called.
- `phase_fcn_generator::F`: callable `(θ, E) -> (phaseE, phaseI)` used to (re)build
    `phase_fcn` whenever the pitch-angle or energy grid changes
- `phase_fcn::P`: tuple `(phaseE, phaseI)` of `[n_θ × n_E]` matrices materialized from the
    generator on the model's scattering θ grid and energy centers.
    Holds 0×0 placeholder matrices until `initialize!(model)` is called.
- `cascading_spec::S`: ionization thresholds and secondary-distribution law used to
    build the cascading transfer matrices
- `cascading_data::C`: cascading transfer matrices (populated lazily by
    `load_or_compute_cascading!`)
"""
mutable struct NeutralSpecies{F, S<:CascadingSpec, C<:SpeciesCascadingCache, P}
    name::Symbol
    density_profile        # untyped: freely replaceable before initialize!(model)
    density::Vector{Float64}
    cross_sections::Matrix{Float64}
    excitation_levels::Matrix{Float64}
    phase_fcn_generator::F
    phase_fcn::P
    cascading_spec::S
    cascading_data::C
end

"""
    NeutralSpecies(name, density_profile; cascading_spec, phase_fcn_generator)

Build a lightweight `NeutralSpecies`. All grid-dependent fields (`density`, `cross_sections`,
`excitation_levels`, `phase_fcn`) start empty and are populated by `initialize!(model)`.
"""
function NeutralSpecies(name::Symbol, density_profile;
                        cascading_spec::CascadingSpec, phase_fcn_generator)
    empty_mat = Matrix{Float64}(undef, 0, 0)
    return NeutralSpecies(
        name,
        density_profile,
        Float64[],
        empty_mat,
        empty_mat,
        phase_fcn_generator,
        (empty_mat, copy(empty_mat)),
        cascading_spec,
        SpeciesCascadingCache(cascading_spec),
    )
end

"""
    apply_density_boundary!(n)

Smoothly drive the top of a density profile to zero so the solver does not see a hard cut.
The last 3 grid points are zeroed and the 3 points below are tapered by an error-function
factor. Negative values (if any sneak in from log-space extrapolation) are clamped to 0.
"""
function apply_density_boundary!(n::AbstractVector)
    n[end-2:end] .= 0
    erf_factor = (erf.((1:-1:-1) / 2) .+ 1) / 2
    n[end-5:end-3] .= erf_factor .* n[end-5:end-3]
    if any(<(0), n)
        @warn "Negative densities found. They were set to 0, but you might want to check " *
              "why that happened"
        n[n .< 0] .= 0
    end
    return n
end

# Per-species variant of load_excitation_threshold (which loads all three species).
function load_excitation_threshold_for(species_name::AbstractString)
    file = pkgdir(AURORA, "internal_data", "data_neutrals", species_name * "_levels.dat")
    return readdlm(file, comments=true, comment_char='%')
end


# ======================================================================================== #
#                        Default species convenience constructors                          #
# ======================================================================================== #

"""
    N2Species(density_profile)
    N2Species(msis_file::AbstractString)

Construct the default N₂ species with cross sections from the built-in `e_N2_*` library,
phase function from [`phase_fcn_N2`](@ref), and cascading described by
[`DefaultCascadingSpecN2`](@ref).

`density_profile` can be an [`MSISDensity`](@ref), a [`VectorDensity`](@ref), or any
callable `h_atm (m) → density (m⁻³)`. Passing an MSIS file path string is a shorthand
for `MSISDensity(msis_file, :N2)`. Grid-dependent fields are populated later by
`initialize!(model)`.
"""
function N2Species(msis_file::AbstractString)
    return N2Species(MSISDensity(msis_file, :N2))
end
function N2Species(density_profile)
    return NeutralSpecies(:N2, density_profile;
                          cascading_spec      = DefaultCascadingSpecN2(),
                          phase_fcn_generator = phase_fcn_N2)
end

"""
    O2Species(density_profile)
    O2Species(msis_file::AbstractString)

Default O₂ species, analogous to [`N2Species`](@ref).
"""
function O2Species(msis_file::AbstractString)
    return O2Species(MSISDensity(msis_file, :O2))
end
function O2Species(density_profile)
    return NeutralSpecies(:O2, density_profile;
                          cascading_spec      = DefaultCascadingSpecO2(),
                          phase_fcn_generator = phase_fcn_O2)
end

"""
    OSpecies(density_profile)
    OSpecies(msis_file::AbstractString)

Default O species, analogous to [`N2Species`](@ref).
"""
function OSpecies(msis_file::AbstractString)
    return OSpecies(MSISDensity(msis_file, :O))
end
function OSpecies(density_profile)
    return NeutralSpecies(:O, density_profile;
                          cascading_spec      = DefaultCascadingSpecO(),
                          phase_fcn_generator = phase_fcn_O)
end

# Read MSIS file, interpolate to grid, return one species' density column.
function load_msis_density(msis_file::AbstractString, species::Symbol, h_atm::AbstractVector)
    msis_raw = load_msis(msis_file)
    msis = interpolate_msis_to_grid(msis_raw.data, h_atm)
    return getproperty(msis, species)
end


# ======================================================================================== #
#                              Display                                                     #
# ======================================================================================== #

function Base.show(io::IO, sp::NeutralSpecies)
    print(io, "NeutralSpecies($(sp.name), ", length(sp.density), " altitudes)")
end

function Base.show(io::IO, ::MIME"text/plain", sp::NeutralSpecies)
    println(io, "NeutralSpecies(", sp.name, "):")
    println(io, "├── Density profile: ", _profile_label(sp.density_profile))
    if isempty(sp.density)
        println(io, "├── (not initialized — call initialize!(model))")
    else
        println(io, "├── Altitudes:        ", length(sp.density))
        println(io, "├── Max density:      ", round(maximum(sp.density), sigdigits=3), " m⁻³")
        println(io, "├── Cross sections:   ", size(sp.cross_sections, 1), " levels × ",
                                              size(sp.cross_sections, 2), " energies")
        println(io, "├── Excitation lvls:  ", size(sp.excitation_levels, 1))
    end
    print(io,   "└── Cascading:        ", sp.cascading_spec.name,
                " (", length(sp.cascading_spec.ionization_thresholds), " thresholds)")
end

_profile_label(d::MSISDensity)  = "MSISDensity($(basename(d.msis_file)), :$(d.species))"
_profile_label(d::VectorDensity) = "VectorDensity($(length(d.h)) points)"
_profile_label(d)               = string(typeof(d))
