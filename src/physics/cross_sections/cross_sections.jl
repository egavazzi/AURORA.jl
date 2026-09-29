"""
    evaluate_in_energy_order(f_sorted, Ep::AbstractVector)

Evaluate a cross-section body `f_sorted`, which assumes energies given in ascending order,
on `Ep` regardless of its order. Returns a floating-point array with the axes of `Ep`, such
that `out[j] == f_sorted([Ep[j]])[1]` for every `j` in `eachindex(Ep)`.
"""
function evaluate_in_energy_order(f_sorted, Ep::AbstractVector)
    # 1-based working copy: the bodies index their input from 1.
    E = Vector{float(eltype(Ep))}(undef, length(Ep))
    copyto!(E, Ep)
    sorted = issorted(E)
    p = sorted ? nothing : sortperm(E)
    σ = f_sorted(sorted ? E : E[p])
    length(σ) == length(E) || throw(DimensionMismatch(
        "cross-section body returned $(length(σ)) values for $(length(E)) energies"))
    sorted || invpermute!(σ, p)
    out = similar(Ep, eltype(σ))
    copyto!(out, σ)
    return out
end

"""
    load_cross_sections(energy_grid)
    load_cross_sections(E_centers::AbstractVector)

Load the cross-sections of the neutrals species for their different energy states.

# Calling
`σ_neutrals = load_cross_sections(energy_grid)`
`σ_neutrals = load_cross_sections(E_centers)`

# Inputs
- `energy_grid`: an `EnergyGrid` object, or
- `E_centers`: energy bin centers (eV). Vector [n\\_E]

# Returns
- `σ_neutrals`: A named tuple containing the cross-sections (m²) for N2, O2, and O.
"""
function load_cross_sections(E_centers::AbstractVector)
    σ_N2 = get_cross_section("N2", E_centers)
    σ_O2 = get_cross_section("O2", E_centers)
    σ_O = get_cross_section("O", E_centers)

    σ_neutrals = (σ_N2 = σ_N2, σ_O2 = σ_O2, σ_O = σ_O)
    return σ_neutrals
end

load_cross_sections(energy_grid::EnergyGrid) = load_cross_sections(energy_grid.E_centers)

"""
    get_cross_section(species, energy_grid)
    get_cross_section(species, E_centers::AbstractVector)

Calculate the cross-section for a given species and their different energy states.

# Calling
`σ_N2 = get_cross_section("N2", energy_grid)`
`σ_N2 = get_cross_section("N2", E_centers)`

# Inputs
- `species`: species name, a `Symbol` or a `String`
- `energy_grid`: an `EnergyGrid` object, or
- `E_centers`: energy bin centers (eV). Vector [n\\_E]

# Returns
- `σ_species`: `[n_levels × n_E]` matrix (m²). Row 1 is the elastic cross section, row
    `i + 1` is `default_channels(species)[i]`.
"""
function get_cross_section(species::Symbol, E_centers::AbstractVector)
    return channel_cross_sections(default_elastic_cross_section(species),
                                  default_channels(species), E_centers)
end

get_cross_section(species::AbstractString, E_centers::AbstractVector) =
    get_cross_section(Symbol(species), E_centers)

get_cross_section(species, energy_grid::EnergyGrid) =
    get_cross_section(species, energy_grid.E_centers)
