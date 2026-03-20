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
