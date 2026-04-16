using DelimitedFiles: readdlm

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

"""
    load_excitation_threshold()

Load the excitation thresholds or energy levels of the different states (vibrational,
rotational, ionization, ...) of the neutrals species as specified in the XX_levels.dat
files. The corresponding names of the states can be found in the XX_levels.name files. XX
refers to N2, O2 or O.

# Calling
`collision_levels = load_excitation_threshold()`

# Returns
- `collision_levels`: A named tuple of matrices, namely `(N2_levels, O2_levels, O_levels)`.
    The matrices have shape [n`_`levels x 2]. The first column contains the energy levels
    and the second column contains the number of secondaries associated to that level (is
    non-zero only for ionized states).
"""
function load_excitation_threshold()
    file_N2_levels = pkgdir(AURORA, "internal_data", "data_neutrals", "N2_levels.dat")
    file_O2_levels = pkgdir(AURORA, "internal_data", "data_neutrals", "O2_levels.dat")
    file_O_levels = pkgdir(AURORA, "internal_data", "data_neutrals", "O_levels.dat")

    N2_levels = readdlm(file_N2_levels, comments=true, comment_char='%')
    O2_levels = readdlm(file_O2_levels, comments=true, comment_char='%')
    O_levels = readdlm(file_O_levels, comments=true, comment_char='%')

    collision_levels = (N2_levels = N2_levels, O2_levels = O2_levels, O_levels = O_levels)
    return collision_levels
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
- `σ_neutrals`: A named tuple containing the cross-sections (m⁻²) for N2, O2, and O.
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
    get_cross_section(species_name, energy_grid)
    get_cross_section(species_name, E_centers::AbstractVector)

Calculate the cross-section for a given species and their different energy states.

# Calling
`σ_N2 = get_cross_section("N2", energy_grid)`
`σ_N2 = get_cross_section("N2", E_centers)`

# Inputs
- `species_name`: name of the species. String
- `energy_grid`: an `EnergyGrid` object, or
- `E_centers`: energy bin centers (eV). Vector [n\\_E]

# Outputs
- `σ_species`: A matrix of cross-section values for each energy state, for the defined
  species
"""
function get_cross_section(species_name, E_centers::AbstractVector)
    state_name = get_level_names(species_name)
    function_name = "e_" * species_name .* state_name

    σ_species = zeros(size(state_name, 1), length(E_centers))
    for i_state in axes(state_name, 1) # loop over the different energy states
        func = getfield(AURORA, Symbol(function_name[i_state])) # get the corresponding function name
        σ_species[i_state, :] .= func(E_centers) # calculate the corresponding cross-section
    end

    return σ_species
end

get_cross_section(species_name, energy_grid::EnergyGrid) = get_cross_section(species_name, energy_grid.E_centers)

"""
    get_level_names(species_name)

Return the names of the excited/ionized states for a given species as a `Vector{String}`.

# Example
```julia
get_level_names("N2")  # → ["_elastic", "_rot0_2", ..., "_ionization"]
```
"""
function get_level_names(species_name)
    filename = pkgdir(AURORA, "internal_data", "data_neutrals", species_name * "_levels.name")
    state_name = readdlm(filename, String, comments=true, comment_char='%')
    return vec(state_name)
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
