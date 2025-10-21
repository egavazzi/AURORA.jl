"""
    interpolate_msis_to_grid(msis_data, h_atm)

Interpolate MSIS data to a custom altitude grid.

This function interpolates all altitude-dependent densities and temperature from MSIS
to a new altitude grid. The interpolation is performed in log space for densities
(which provides exponential extrapolation) and linear space for temperature.

# Arguments
- `msis_data`: NamedTuple from `load_msis_data()` or `load_msis().data`
- `h_atm`: Target altitude grid (m). Vector [nZ]

# Returns
A NamedTuple with the following fields (all vectors of length nZ):
- `height_km`: converted target altitude grid (km)
- `N2`: interpolated N₂ density (m⁻³)
- `O2`: interpolated O₂ density (m⁻³)
- `O`: interpolated O density (m⁻³)
- `T`: interpolated temperature (K)

# Notes
- Densities are interpolated in log space for better exponential behavior
- Temperature is interpolated in linear space
"""
function interpolate_msis_to_grid(msis_data, h_atm)
    # Extract altitude from MSIS data
    z_msis = msis_data.height_km

    # Convert target altitude to km
    h_atm_km = h_atm / 1e3

    # Interpolate densities (log space)
    N2 = interpolate_profile(msis_data.N2, z_msis, h_atm; log_interpolation = true)
    O2 = interpolate_profile(msis_data.O2, z_msis, h_atm; log_interpolation = true)
    O = interpolate_profile(msis_data.O, z_msis, h_atm; log_interpolation = true)
    He = interpolate_profile(msis_data.He, z_msis, h_atm; log_interpolation = true)
    H = interpolate_profile(msis_data.H, z_msis, h_atm; log_interpolation = true)
    Ar = interpolate_profile(msis_data.Ar, z_msis, h_atm; log_interpolation = true)
    N = interpolate_profile(msis_data.N, z_msis, h_atm; log_interpolation = true)
    anomalousO = interpolate_profile(msis_data.anomalousO, z_msis, h_atm; log_interpolation = true)
    NO = interpolate_profile(msis_data.NO, z_msis, h_atm; log_interpolation = true)

    # Interpolate temperature (linear space)
    T = interpolate_profile(msis_data.T, z_msis, h_atm; log_interpolation = false)

    return (; height_km = h_atm_km, N2, O2, O, He, H, Ar, N, anomalousO, NO, T)
end
