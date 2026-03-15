using DataInterpolations: PCHIPInterpolation, ExtrapolationType
using DelimitedFiles: readdlm
using SpecialFunctions: erf


"""
    setup(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith=0)

Load the atmosphere, the energy grid, and the collision data.

## Calling
`state = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)`

## Inputs
- `altitude_lims`: the altitude limits, in km, for the bottom and top of the ionosphere in our simulation
- `θ_lims`: pitch-angle limits of the electron beams (e.g. 180:-10:0), where 180°
    corresponds to field aligned down, and 0° field aligned up. Vector [n_beam]
- `E_max`: upper limit for the energy grid (in eV)
- `msis_file`: path to the msis file to use
- `iri_file`: path to the iri file to use
- `B_angle_to_zenith`: angle between magnetic field and zenith (degrees)

## Returns
A NamedTuple with fields:
- `altitude_grid::AltitudeGrid`
- `energy_grid::EnergyGrid`
- `pitch_angle_grid::PitchAngleGrid`
- `scattering::ScatteringData`
- `ionosphere::Ionosphere`
- `cross_sections::CrossSectionData`
- `B_angle_to_zenith::Real`
- `s_field::Vector`
"""
function setup(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith=0)
    altitude_grid = AltitudeGrid(altitude_lims[1], altitude_lims[2])
    energy_grid = EnergyGrid(E_max)
    pitch_angle_grid = PitchAngleGrid(θ_lims)
    scattering = ScatteringData(θ_lims)
    ionosphere = Ionosphere(msis_file, iri_file, altitude_grid.h)
    cross_sections = CrossSectionData(energy_grid)
    s_field = altitude_grid.h ./ cosd(B_angle_to_zenith)

    return (; altitude_grid, energy_grid, pitch_angle_grid,
              scattering, ionosphere, cross_sections,
              B_angle_to_zenith, s_field)
end



"""
    make_altitude_grid(bottom_altitude, top_altitude; dz_max=25)

Create an altitude grid based on the altitude limits given as input. It uses
constant steps of 150m for altitudes below 100km, and a non-linear grid above 100km.

# Calling
`h_atm = make_altitude_grid(bottom_altitude, top_altitude; max_dz=25)`

# Arguments
- `bottom_altitude`: altitude, in km, for the bottom of the simulation
- `top_altitude`: altitude, in km, for the top of the simulation

# Keyword Arguments
- `dz_max = 25`: maximum step size, in km. Relevant for high altitudes where the numbers
    can get large.

# Outputs
- `h_atm`: altitude (m). Vector [nZ]
"""
function make_altitude_grid(bottom_altitude, top_altitude; dz_max = 25)
    Δz(n) = 150 .+
            150 / 200 * (0:(n - 1)) .+
            1.2 * exp.(Complex.(((0:(n - 1)) .- 150) / 22) .^ 0.9)
    dz = real.(Δz(500))
    dz = min.(dz, dz_max * 1e3)
    # With the default dz_max of 25 km, the grid can go up to ~2700 km (which is way
    # too high, don't do that).
    h_atm = 100e3 .+ cumsum(dz) .- dz[1]
    h_atm = vcat((bottom_altitude * 1e3):(h_atm[2] - h_atm[1]):h_atm[1], h_atm[2:end]) # add altitude steps under 100km
    i_zmax = findlast(h_atm .<= top_altitude * 1e3)
    h_atm = h_atm[1:i_zmax]
    return h_atm
end



"""
    make_energy_grid(E_max)

Create an energy grid based on the maximum energy `E_max` given as input.

# Calling
`E_edges, E_centers, ΔE = make_energy_grid(E_max)`

# Inputs
- `E_max`: upper limit for the energy grid (in eV)

# Outputs
- `E_edges`: energy bin edges (eV). Vector [nE + 1] (includes the last upper edge)
- `E_centers`: energy bin centers (eV). Vector [nE]
- `ΔE`: energy bin widths (eV). Vector [nE]
"""
function make_energy_grid(E_max)
    if E_max > 1e6
        error("AURORA does not support energies above 1 MeV. Please use a lower value for E_max.")
    end
    E_function(X, dE_initial, dE_final, C, X0) = dE_initial + (1 + tanh(C * (X - X0))) / 2 * dE_final
    E = cumsum(E_function.(0:100000, 0.15, 11.5, 0.05, 80)) .+ 1.9
    iE_max = findmin(abs.(E .- E_max))[2];  # find the index for the upper limit of the energy grid
    E = E[1:iE_max];                        # crop E accordingly
    ΔE = diff(E); ΔE = [ΔE; ΔE[end]]
    E_edges = [E; E[end] + ΔE[end]]         # add the last upper edge
    E_centers = E_edges[1:end-1] .+ ΔE ./ 2
    return E_edges, E_centers, ΔE
end



"""
    load_scattering_matrices(θ_lims)

Load the scattering matrices for the given pitch-angle limits.

# Calling
`μ_lims, μ_center, scattering = load_scattering_matrices(θ_lims)`

# Inputs
- `θ_lims`: pitch angle limits of the e- beams (deg). Vector [n_beam + 1]

# Outputs
- `μ_lims`: cosine of the pitch angle limits of the e- beams. Vector [n_beam + 1]
- `μ_center`: cosine of the pitch angle of the middle of the e- beams. Vector [n_beam]
- `scattering`: Tuple with several of the scattering informations, namely scattering
    = `(P_scatter, Ω_beam_relative, Ω_beam)`
    + `P_scatter`: probabilities for scattering in 3D from beam to beam. Matrix [n`_`direction x n`_`direction]
    + `Ω_beam_relative`: relative contribution from within each beam. Matrix [n`_`beam x n`_`direction]
    + `Ω_beam`: solid angle for each stream (ster). Vector [n_beam]
    + `θ_scatter`: scattering angles used in the calculations. Vector [n_direction]
"""
function load_scattering_matrices(θ_lims)
    validate_θ_lims(θ_lims)
    μ_lims = cosd.(θ_lims);
    μ_center = mu_avg(θ_lims);
    Ω_beam = beam_weight(θ_lims); # this beam weight is calculated in a continuous way
    P_scatter, _, Ω_beam_relative, θ₁ = find_scattering_matrices(θ_lims, 720)
    scattering = (P_scatter = P_scatter, Ω_beam_relative = Ω_beam_relative,
                     Ω_beam = Ω_beam, θ_scatter = θ₁)

    return μ_lims, μ_center, scattering
end



"""
    validate_θ_lims(θ_lims)

Validate that the pitch-angle limits `θ_lims` are correctly specified.
Throws an `ArgumentError` if:
- `θ_lims` does not include 180° (field-aligned downward)
- `θ_lims` does not include 0° (field-aligned upward)
- `θ_lims` is not in descending order

# Inputs
- `θ_lims`: pitch-angle limits of the electron beams (e.g. 180:-10:0)
"""
function validate_θ_lims(θ_lims)
    if maximum(θ_lims) != 180
        throw(ArgumentError("θ_lims must include 180° (field-aligned downward). " *
              "Got maximum of $(maximum(θ_lims))°. " *
              "Example of valid input: 180:-10:0."))
    end
    if minimum(θ_lims) != 0
        throw(ArgumentError("θ_lims must include 0° (field-aligned upward). " *
              "Got minimum of $(minimum(θ_lims))°. " *
              "Example of valid input: 180:-10:0."))
    end
    if !issorted(θ_lims, rev=true)
        throw(ArgumentError("θ_lims must be in descending order (e.g., 180:-10:0). " *
              "Got: $θ_lims"))
    end
    return nothing
end



"""
    load_neutral_densities(msis_file, h_atm)

Load the neutral densities and temperature from a MSIS file that was generated and saved
using AURORA's MSIS interface. Then interpolate the profiles over AURORA's altitude grid.

Upper boundary conditions are applied to smoothly transition the densities to zero.

# Calling
`n_neutrals, Tn = load_neutral_densities(msis_file, h_atm)`

# Inputs
- `msis_file`: absolute path to the msis file to read n_neutrals and Tn from. String
- `h_atm`: altitude (m). Vector [nZ]

# Returns
- `n_neutrals`: neutral densities (m⁻³). Named tuple of vectors ([nZ], ..., [nZ])
- `Tn`: neutral temperature (K). Vector [nZ]

# See also
[`load_msis`](@ref), [`interpolate_msis_to_grid`](@ref)
"""
function load_neutral_densities(msis_file, h_atm)
    msis_raw = load_msis(msis_file)
    msis = interpolate_msis_to_grid(msis_raw.data, h_atm)
    nN2 = msis.N2
    nO2 = msis.O2
    nO = msis.O

    # Apply upper boundary conditions to smoothly transition densities to zero
    # Set the last 3 points to zero
    nN2[end-2:end] .= 0
    nO2[end-2:end] .= 0
    nO[end-2:end] .= 0

    # Apply smooth transition using error function for the 3 points below
    erf_factor = (erf.((1:-1:-1) / 2) .+ 1) / 2
    nN2[end-5:end-3] .= erf_factor .* nN2[end-5:end-3]
    nO2[end-5:end-3] .= erf_factor .* nO2[end-5:end-3]
    nO[end-5:end-3] .= erf_factor .* nO[end-5:end-3]

    # Ensure no negative densities
    if any(nN2 .< 0) || any(nO2 .< 0) || any(nO .< 0)
        @warn "Negative densities found. They were set to 0, but you might want to check " *
              "why that happened"
        nN2[nN2 .< 0] .= 0
        nO2[nO2 .< 0] .= 0
        nO[nO .< 0] .= 0
    end

    n_neutrals = (; nN2, nO2, nO)
    return n_neutrals, msis.T
end



"""
    load_electron_densities(iri_file, h_atm)

Load the electron density and temperature from an IRI file that was generated and saved
using AURORA's IRI interface. Then interpolate the profiles over AURORA's altitude grid.

# Calling
`ne, Te = load_electron_densities(iri_file, h_atm)`

# Inputs
- `iri_file`: absolute path to the iri file to read ne and Te from. String
- `h_atm`: altitude (m). Vector [nZ]

# Returns
- `ne`: e- density (m⁻³). Vector [nZ]
- `Te`: e- temperature (K). Vector [nZ]

# See also
[`load_iri`](@ref), [`interpolate_iri_to_grid`](@ref)
"""
function load_electron_densities(iri_file, h_atm)
    iri_raw = load_iri(iri_file)
    iri = interpolate_iri_to_grid(iri_raw.data, h_atm)

    return iri.ne, iri.Te
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
    filename =  pkgdir(AURORA, "internal_data", "data_neutrals", species_name * "_levels.name")
    state_name = readdlm(filename, String, comments=true, comment_char='%')
    function_name = "e_" * species_name .* state_name

    σ_species = zeros(size(state_name, 1), length(E_centers))
    for i_state in axes(state_name, 1) # loop over the different energy states
        func = getfield(AURORA, Symbol(function_name[i_state])) # get the corresponding function name
        σ_species[i_state, :] .= func(E_centers) # calculate the corresponding cross-section
    end

    return σ_species
end

get_cross_section(species_name, energy_grid::EnergyGrid) = get_cross_section(species_name, energy_grid.E_centers)
