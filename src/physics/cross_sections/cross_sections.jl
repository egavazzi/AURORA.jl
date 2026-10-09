"""
    zero_below_energy_loss!(σ, E_centers, E_loss, label)

Set `σ` to zero in the energy bins whose center is below `E_loss`, where the collision
cannot happen, and warn (once per session for each `label`) when any of those values was
non-zero: it points to a cross-section fit that starts below its energy loss.
"""
function zero_below_energy_loss!(σ, E_centers, E_loss, label)
    below = (E_centers .< E_loss) .& (σ .> 0)
    if any(below)
        @warn "$(label) has a non-zero cross section below its energy loss of $(E_loss) eV, \
               at bin centers $(round.(E_centers[below]; digits = 3)) eV. These values are set \
               to zero." maxlog = 1 _id = (label, E_loss)
        σ[below] .= 0
    end
    return σ
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
                                  default_channels(species), E_centers;
                                  species_name = String(species))
end

get_cross_section(species::AbstractString, E_centers::AbstractVector) =
    get_cross_section(Symbol(species), E_centers)

get_cross_section(species, energy_grid::EnergyGrid) =
    get_cross_section(species, energy_grid.E_centers)
