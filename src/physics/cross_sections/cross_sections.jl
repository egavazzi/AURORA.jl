
"""
    evaluate_in_energy_order(f_sorted, Ep::AbstractVector)

Evaluate a cross-section body `f_sorted`, which assumes energies given in ascending order,
on `Ep` regardless of its order. Returns a floating-point array with the axes of `Ep`, such
that `out[j] == f_sorted([Ep[j]])[1]` for every `j` in `eachindex(Ep)`.
"""
function evaluate_in_energy_order(f_sorted, Ep::AbstractVector)
    idxs = collect(eachindex(Ep))
    E = collect(float(eltype(Ep)), Ep) # E[k] is the k-th element of Ep, i.e. Ep[idxs[k]]

    if issorted(E)
        σ = f_sorted(E)
        out = similar(Ep, eltype(σ))
        for k in eachindex(idxs)
            out[idxs[k]] = σ[k]
        end
        return out
    end

    p = sortperm(E)
    σ_sorted = f_sorted(E[p])
    out = similar(Ep, eltype(σ_sorted))
    for k in eachindex(p)
        out[idxs[p[k]]] = σ_sorted[k]
    end
    return out
end

"""
    load_cross_sections(energy_grid)
    load_cross_sections(E_centers::AbstractVector)

Cross sections of the built-in N₂, O₂ and O channel tables, evaluated on the energy grid.

# Inputs
- `energy_grid`: an `EnergyGrid` object, or
- `E_centers`: energy bin centers (eV). Vector [n\\_E]

# Returns
- `σ_neutrals`: named tuple `(σ_N2, σ_O2, σ_O)` of `[n_levels × n_E]` matrices (m²).
"""
function load_cross_sections(E_centers::AbstractVector)
    return (σ_N2 = get_cross_section(:N2, E_centers),
            σ_O2 = get_cross_section(:O2, E_centers),
            σ_O  = get_cross_section(:O,  E_centers))
end

load_cross_sections(energy_grid::EnergyGrid) = load_cross_sections(energy_grid.E_centers)

"""
    get_cross_section(species, energy_grid)
    get_cross_section(species, E_centers::AbstractVector)

Cross sections of one of the built-in species (`:N2`, `:O2` or `:O`), evaluated on the energy
grid.

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
