using DataInterpolations: PCHIPInterpolation, ExtrapolationType
using DelimitedFiles: readdlm
using SpecialFunctions: erf

"""
    setup(altitude_lims, θ_lims, E_max, msis_file, iri_file)

Load the atmosphere, the energy grid, the collision cross-sections, ...

## Calling
`h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, θ_lims, μ_lims,
μ_center, μ_scatterings = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file)`

## Inputs
- `altitude_lims`: the altitude limits, in km, for the bottom and top of the ionosphere in our simulation
- `θ_lims`: pitch-angle limits of the electron beams (e.g. 180:-10:0), where 180°
    corresponds to field aligned down, and 0° field aligned up. Vector [n_beam]
- `E_max`: upper limit for the energy grid (in eV)
- `msis_file`: path to the msis file to use
- `iri_file`: path to the iri file to use

## Outputs
- `h_atm`: altitude (m). Vector [nZ]
- `ne`: e- density (m⁻³). Vector [nZ]
- `Te`: e- temperature (K). Vector [nZ]
- `Tn`: neutral temperature (K). Vector [nZ]
- `E`: energy grid (eV). Vector [nE]
- `dE`: energy bin sizes(eV). Vector [nE]
- `n_neutrals`: neutral densities (m⁻³). Named tuple of vectors ([nZ], ..., [nZ])
- `E_levels_neutrals`: collisions energy levels and number of secondary e- produced. Named
    tuple of matrices ([n`_`levels x 2], ..., [n`_`levels x 2])
- `σ_neutrals`: collision cross-sections (m⁻²). Named tuple of matrices ([n`_`levels x nE],
    ..., [n`_`levels x nE])
- `θ_lims`: pitch angle limits of the e- beams (deg). Vector [n_beam + 1]
- `μ_lims`: cosine of the pitch angle limits of the e- beams. Vector [n_beam + 1]
- `μ_center`: cosine of the pitch angle of the middle of the e- beams. Vector [n_beam]
- `μ_scatterings`: Named tuple with several of the scattering informations, namely
    μ`_`scatterings = `(Pmu2mup, BeamWeight_relative, BeamWeight)`
    + `Pmu2mup`: probabilities for scattering in 3D from beam to beam. Matrix [n`_`direction x n`_`direction]
    + `BeamWeight_relative`: relative contribution from within each beam. Matrix [n`_`beam x n`_`direction]
    + `BeamWeight`: solid angle for each stream (ster). Vector [n_beam]
    + `theta1`: scattering angles used in the calculations. Vector [n_direction]
"""
function setup(altitude_lims, θ_lims, E_max, msis_file, iri_file)
    h_atm = make_altitude_grid(altitude_lims[1], altitude_lims[2])
    E, dE = make_energy_grid(E_max)
    μ_lims, μ_center, μ_scatterings = load_scattering_matrices(θ_lims)
    n_neutrals, Tn = load_neutral_densities(msis_file, h_atm)
    ne, Te = load_electron_densities(iri_file, h_atm)
    E_levels_neutrals = load_excitation_threshold()
    σ_neutrals = load_cross_sections(E, dE)

    return h_atm, ne, Te, Tn, E, dE,
    n_neutrals, E_levels_neutrals, σ_neutrals,
    μ_lims, μ_center, μ_scatterings
end



"""
    make_altitude_grid(bottom_altitude, top_altitude)

Create an altitude grid based on the altitude limits given as input. It uses
constant steps of 150m for altitudes below 100km, and a non-linear grid above 100km.

# Calling
`h_atm = make_altitude_grid(bottom_altitude, top_altitude)

# Inputs
- `top_altitude`: the altitude, in km, for the top of the ionosphere in our simulation
- `bottom_altitude`: the altitude, in km, for the bottom of the ionosphere in our simulation

# Outputs
- `h_atm`: altitude (m). Vector [nZ]
"""
function make_altitude_grid(bottom_altitude, top_altitude)
    Δz(n) = 150 .+
            150 / 200 * (0:(n - 1)) .+
            1.2 * exp.(Complex.(((0:(n - 1)) .- 150) / 22) .^ .9)
    h_atm = 100e3 .+ cumsum(real.(Δz(500))) .- real.(Δz(1))
    h_atm = vcat((bottom_altitude * 1e3):(h_atm[2] - h_atm[1]):h_atm[1], h_atm[2:end]) # add altitude steps under 100km
    i_zmax = findmin(abs.(h_atm .- top_altitude * 1e3))[2]
    h_atm = h_atm[1:i_zmax]
    return h_atm
end



"""
    make_energy_grid(E_max)

Create an energy grid based on the maximum energy `E_max` given as input.

# Calling
`E, dE = make_energy_grid(E_max)`

# Inputs
- `E_max`: upper limit for the energy grid (in eV)

# Outputs
- `E`: energy grid (eV). Vector [nE]
- `dE`: energy bin sizes(eV). Vector [nE]
"""
function make_energy_grid(E_max)
    E_function(X, dE_initial, dE_final, C, X0) = dE_initial + (1 + tanh(C * (X - X0))) / 2 * dE_final
    E = cumsum(E_function.(0:4000, 0.15, 11.5, 0.05, 80)) .+ 1.9
    iE_max = findmin(abs.(E .- E_max))[2];  # find the index for the upper limit of the energy grid
    E = E[1:iE_max];                        # crop E accordingly
    dE = diff(E); dE = [dE; dE[end]]
    return E, dE
end



"""
    load_scattering_matrices(θ_lims)

Create an energy grid based on the maximum energy `E_max` given as input.

# Calling
`μ_lims, μ_center, μ_scatterings = load_scattering_matrices(θ_lims)`

# Inputs
- `θ_lims`: pitch angle limits of the e- beams (deg). Vector [n_beam + 1]

# Outputs
- `μ_lims`: cosine of the pitch angle limits of the e- beams. Vector [n_beam + 1]
- `μ_center`: cosine of the pitch angle of the middle of the e- beams. Vector [n_beam]
- `μ_scatterings`: Tuple with several of the scattering informations, namely μ`_`scatterings
    = `(Pmu2mup, BeamWeight_relative, BeamWeight)`
    + `Pmu2mup`: probabilities for scattering in 3D from beam to beam. Matrix [n`_`direction x n`_`direction]
    + `BeamWeight_relative`: relative contribution from within each beam. Matrix [n`_`beam x n`_`direction]
    + `BeamWeight`: solid angle for each stream (ster). Vector [n_beam]
    + `theta1`: scattering angles used in the calculations. Vector [n_direction]
"""
function load_scattering_matrices(θ_lims)
    μ_lims = cosd.(θ_lims);
    μ_center = mu_avg(θ_lims);
    BeamWeight = beam_weight(θ_lims); # this beam weight is calculated in a continuous way
    Pmu2mup, _, BeamWeight_relative, θ₁ = find_scattering_matrices(θ_lims, 720)
    μ_scatterings = (Pmu2mup = Pmu2mup, BeamWeight_relative = BeamWeight_relative,
                     BeamWeight = BeamWeight, theta1 = θ₁)

    return μ_lims, μ_center, μ_scatterings
end



"""
    load_neutral_densities(msis_file, h_atm)

Load the neutral densities and temperature.

# Calling
`n_neutrals, Tn = load_neutral_densities(msis_file, h_atm)`

# Inputs
- `msis_file`: absolute path to the msis file to read n_neutrals and Tn from. String
- `h_atm`: altitude (m). Vector [nZ]

# Returns
- `n_neutrals`: neutral densities (m⁻³). Named tuple of vectors ([nZ], ..., [nZ])
- `Tn`: neutral temperature (K). Vector [nZ]
"""
function load_neutral_densities(msis_file, h_atm)
    data_msis = readdlm(msis_file, skipstart=14) # read the file without the headers
    data_msis[isnan.(data_msis)] .= 0

    z_msis = data_msis[:, 1] # extract the z-grid of the msis data
    data_msis = data_msis[:, [1:5; end]] # keep only data of interest (and also avoid NaN values)

    # create the interpolator
    msis_interpolator = [PCHIPInterpolation(data_msis[:, i], z_msis;
                                            extrapolation = ExtrapolationType.Extension)
                         for i in axes(data_msis, 2)]
    # interpolate the msis data over our h_atm grid
    msis_interpolated = [msis_interpolator[i](h_atm / 1e3) for i in axes(data_msis, 2)]

    # extract the neutral densities
    nN2 = msis_interpolated[3] # already in m⁻³
    nO2 = msis_interpolated[4] # already in m⁻³
    nO = msis_interpolated[5]  # already in m⁻³
    Tn = msis_interpolated[6]

    nO[end-2:end] .= 0
	nN2[end-2:end] .= 0
	nO2[end-2:end] .= 0
    erf_factor = (erf.((1:-1:-1) / 2) .+ 1) / 2
	nO[end-5:end-3] .= erf_factor .* nO[end-5:end-3]
	nN2[end-5:end-3] .= erf_factor .* nN2[end-5:end-3]
	nO2[end-5:end-3] .= erf_factor .* nO2[end-5:end-3]

    n_neutrals = (nN2 = nN2, nO2 = nO2, nO = nO)
    return n_neutrals, Tn
end



"""
    load_electron_densities(iri_file, h_atm)

Load the electron density and temperature.

# Calling
`ne, Te = load_electron_densities(iri_file, h_atm)`

# Inputs
- `iri_file`: absolute path to the iri file to read ne and Te from. String
- `h_atm`: altitude (m). Vector [nZ]

# Returns
- `ne`: e- density (m⁻³). Vector [nZ]
- `Te`: e- temperature (K). Vector [nZ]
"""
function load_electron_densities(iri_file, h_atm)
    # read the file and extract z-grid of the iri data
    data_iri = readdlm(iri_file, skipstart=14)
    data_iri[isnan.(data_iri)] .= 0

    z_iri = data_iri[:, 1]

    # create the interpolator
    iri_interpolator = [PCHIPInterpolation(data_iri[:, i], z_iri;
                                           extrapolation = ExtrapolationType.Extension)
                        for i in axes(data_iri, 2)]
    # interpolate the iri data over our h_atm grid
    iri_interpolated = [iri_interpolator[i](h_atm / 1e3) for i in axes(data_iri, 2)]

    # extract electron density and temperature
    ne = iri_interpolated[2] # from cm⁻³ to m⁻³
    Te = iri_interpolated[5]

    ne[ne .< 0] .= 1
    Te[Te .== -1] .= 350

    return ne, Te
end



"""
    load_excitation_threshold()

Load the excitation thresholds or energy levels of the different states (vibrational,
rotational, ionization, ...) of the neutrals species as specified in the XX_levels.dat
files. The corresponding names of the states can be found in the XX_levels.name files. XX
refers to N2, O2 or O.

# Calling
`E_levels_neutrals = load_excitation_threshold()`

# Returns
- `E_levels_neutrals`: A named tuple of matrices, namely `(N2_levels, O2_levels, O_levels)`.
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

	E_levels_neutrals = (N2_levels = N2_levels, O2_levels = O2_levels, O_levels = O_levels)
    return E_levels_neutrals
end



"""
    load_cross_sections(E, dE)

Load the cross-sections of the neutrals species for their different energy states.

# Calling
`σ_neutrals = load_cross_sections(E, dE)`

# Inputs
- `E`: energy grid (eV). Vector [nE]
- `dE`: energy grid step size (eV). Vector [nE]

# Returns
- `σ_neutrals`: A named tuple containing the cross-sections (m⁻²) for N2, O2, and O.
"""
function load_cross_sections(E, dE)
    σ_N2 = get_cross_section("N2", E, dE)
    σ_O2 = get_cross_section("O2", E, dE)
    σ_O = get_cross_section("O", E, dE)

    σ_neutrals = (σ_N2 = σ_N2, σ_O2 = σ_O2, σ_O = σ_O)
    return σ_neutrals
end



"""
    get_cross_section(species_name, E, dE)

Calculate the cross-section for a given species and their different energy states.

# Calling
`σ_N2 = get_cross_section("N2", E, dE)`

# Inputs
- `species_name`: name of the species. String
- `E`: energy grid (eV). Vector [nE]
- `dE`: energy grid step size (eV). Vector [nE]

# Outputs
- `σ_species`: A matrix of cross-section values for each energy state, for the defined
  species
"""
function get_cross_section(species_name, E, dE)
    filename =  pkgdir(AURORA, "internal_data", "data_neutrals", species_name * "_levels.name")
    state_name = readdlm(filename, String, comments=true, comment_char='%')
    function_name = "e_" * species_name .* state_name

    σ_species = zeros(size(state_name, 1), length(E))
    for i_state in axes(state_name, 1) # loop over the different energy states
        func = getfield(AURORA, Symbol(function_name[i_state])) # get the corresponding function name
        σ_species[i_state, :] .= func(E .+ dE/2) # calculate the corresponding cross-section
        # Note that we use `E .+ dE/2` as input to the cross-section functions, as this
        # corresponds to the middle energy of the energy bins (E is the energy *grid*).
    end

    return σ_species
end
