using DataInterpolations: PCHIPInterpolation, ExtrapolationType

"""
    interpolate_iri_to_grid(iri_data, h_atm)

Interpolate IRI data to a custom altitude grid.

This function interpolates all altitude-dependent densities and temperatures from IRI
to a new altitude grid, while preserving scalar parameters (peak densities/heights, TEC, etc.).

# Arguments
- `iri_data`: NamedTuple from `load_iri_data()` or `load_iri().data`
- `h_atm`: Target altitude grid (m). Vector [nZ]

# Returns
A NamedTuple with the same fields as `iri_data`, but with interpolated profiles:
- `height_km`: converted target altitude grid (km). Vector [nZ]
- `ne`: interpolated electron density (m⁻³). Vector [nZ]
- `Tn`: interpolated neutral temperature (K). Vector [nZ]
- `Ti`: interpolated ion temperature (K). Vector [nZ]
- `Te`: interpolated electron temperature (K). Vector [nZ]
- `nO⁺`, `nH⁺`, `nHe⁺`, etc.: interpolated ion densities (m⁻³). Vector [nZ]
- `NmF2`, `hmF2`, `NmF1`, etc.: preserved scalar parameters (unchanged from input)
"""
function interpolate_iri_to_grid(iri_data, h_atm)
    # Extract altitude from IRI data
    z_iri = iri_data.height_km

    # Convert target altitude to km for output
    h_atm_km = h_atm / 1e3

    # Interpolate electron density (log space)
    ne = interpolate_profile(iri_data.ne, z_iri, h_atm; log_interpolation = true)

    # Interpolate temperatures (linear space)
    Tn = interpolate_profile(iri_data.Tn, z_iri, h_atm; log_interpolation = false)
    Ti = interpolate_profile(iri_data.Ti, z_iri, h_atm; log_interpolation = false)
    Te = interpolate_profile(iri_data.Te, z_iri, h_atm; log_interpolation = false)

    # Interpolate ion densities (log space)
    nO⁺ = interpolate_profile(iri_data.nO⁺, z_iri, h_atm; log_interpolation = true)
    nH⁺ = interpolate_profile(iri_data.nH⁺, z_iri, h_atm; log_interpolation = true)
    nHe⁺ = interpolate_profile(iri_data.nHe⁺, z_iri, h_atm; log_interpolation = true)
    nO2⁺ = interpolate_profile(iri_data.nO2⁺, z_iri, h_atm; log_interpolation = true)
    nNO⁺ = interpolate_profile(iri_data.nNO⁺, z_iri, h_atm; log_interpolation = true)
    nN⁺ = interpolate_profile(iri_data.nN⁺, z_iri, h_atm; log_interpolation = true)

    # Handle cluster ions specially: IRI-2016 outputs negative values for some reason (it
    # seens to be hard-coded in the source code to be equal to -1% the electron density??)
    # Interpolate absolute values to get positive densities.
    nCI = interpolate_profile(abs.(iri_data.nCI), z_iri, h_atm; log_interpolation = true)

    # Return new named tuple with interpolated profiles and preserved scalar parameters
    return (height_km = h_atm_km,
            ne = ne,
            Tn = Tn,
            Ti = Ti,
            Te = Te,
            nO⁺ = nO⁺,
            nH⁺ = nH⁺,
            nHe⁺ = nHe⁺,
            nO2⁺ = nO2⁺,
            nNO⁺ = nNO⁺,
            nCI = nCI,
            nN⁺ = nN⁺,
            NmF2 = iri_data.NmF2,
            hmF2 = iri_data.hmF2,
            NmF1 = iri_data.NmF1,
            hmF1 = iri_data.hmF1,
            NmE = iri_data.NmE,
            hmE = iri_data.hmE,
            TEC = iri_data.TEC,
            EqVertIonDrift = iri_data.EqVertIonDrift,
            foF2 = iri_data.foF2)
end

"""
    interpolate_profile(data_values, data_altitude_km, target_altitude_m;
                       log_interpolation=true)

Core interpolation function for atmospheric profiles (densities, temperatures, etc.).

This function interpolates data from one altitude grid to another using PCHIP interpolation.
It can handle both linear and logarithmic interpolation.

# Arguments
- `data_values::Vector`: The data to interpolate (e.g., density or temperature)
- `data_altitude_km::Vector`: Altitude grid of the input data (km)
- `target_altitude_m::Vector`: Target altitude grid for interpolation (m)

# Keyword Arguments
- `log_interpolation::Bool=true`: If `true`, interpolation is done in log space (exponential
    extrapolation). Recommended for densities. Use `false` for temperatures.

# Returns
- `interpolated::Vector`: Interpolated data on the target altitude grid
"""
function interpolate_profile(data_values, data_altitude_km, target_altitude_m;
                             log_interpolation = true)
    # Convert target altitude from meters to kilometers
    target_altitude_km = target_altitude_m / 1e3

    # Prepare data for interpolation
    if log_interpolation
        # Use log interpolation for exponential-like profiles (densities)
        interpolator = PCHIPInterpolation(log.(data_values), data_altitude_km;
                                          extrapolation = ExtrapolationType.Linear)
        interpolated = exp.(interpolator(target_altitude_km))
    else
        # Use linear interpolation (e.g., for temperatures)
        interpolator = PCHIPInterpolation(data_values, data_altitude_km;
                                          extrapolation = ExtrapolationType.Linear)
        interpolated = interpolator(target_altitude_km)
    end

    # Check for negative values and throw error if found
    if any(interpolated .< 0)
        error("Interpolation resulted in negative values. Check input data or interpolation method.")
    end

    return interpolated
end
