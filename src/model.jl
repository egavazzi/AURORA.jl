"""
    AuroraModel{AG, EG, PAG, SD, IO, CS, FT, V}

Container for the grids, atmosphere, and collision data used by an AURORA simulation.
"""
struct AuroraModel{AG<:AltitudeGrid, EG<:EnergyGrid, PAG<:PitchAngleGrid,
                   SD<:ScatteringData, IO<:Ionosphere, CS<:CrossSectionData,
                   FT, V<:AbstractVector{FT}}
    altitude_grid::AG
    energy_grid::EG
    pitch_angle_grid::PAG
    scattering::SD
    ionosphere::IO
    cross_sections::CS
    B_angle_to_zenith::FT
    s_field::V
end

"""
    AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith=0;
                verbose=true)

Load the atmosphere, the energy grid, and the collision data into an `AuroraModel`.

## Calling
`model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)`

## Inputs
- `altitude_lims`: the altitude limits, in km, for the bottom and top of the ionosphere in our simulation
- `θ_lims`: pitch-angle limits of the electron beams (e.g. 180:-10:0), where 180°
    corresponds to field aligned down, and 0° field aligned up. Vector [n_beam]
- `E_max`: upper limit for the energy grid (in eV)
- `msis_file`: path to the msis file to use
- `iri_file`: path to the iri file to use
- `B_angle_to_zenith`: angle between magnetic field and zenith (degrees)
- `verbose=true`: print progress messages when loading stuff

## Returns
An `AuroraModel` with fields:
- `altitude_grid::AltitudeGrid`
- `energy_grid::EnergyGrid`
- `pitch_angle_grid::PitchAngleGrid`
- `scattering::ScatteringData`
- `ionosphere::Ionosphere`
- `cross_sections::CrossSectionData`
- `B_angle_to_zenith::Real`
- `s_field::Vector`
"""
function AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith=0;
                     verbose=true)
    altitude_grid = AltitudeGrid(altitude_lims[1], altitude_lims[2])
    energy_grid = EnergyGrid(E_max)
    pitch_angle_grid = PitchAngleGrid(θ_lims)
    scattering = ScatteringData(pitch_angle_grid; verbose)
    ionosphere = Ionosphere(msis_file, iri_file, altitude_grid.h)
    cross_sections = CrossSectionData(energy_grid)
    s_field = altitude_grid.h ./ cosd(B_angle_to_zenith)

    FT = promote_type(eltype(s_field), typeof(B_angle_to_zenith))

    return AuroraModel{typeof(altitude_grid), typeof(energy_grid), typeof(pitch_angle_grid),
                       typeof(scattering), typeof(ionosphere), typeof(cross_sections),
                       FT, typeof(s_field)}(
        altitude_grid, energy_grid, pitch_angle_grid,
        scattering, ionosphere, cross_sections,
        FT(B_angle_to_zenith), s_field
    )
end

function Base.show(io::IO, model::AuroraModel)
    print(io, "AuroraModel($(model.altitude_grid), $(model.energy_grid))")
end

function Base.show(io::IO, ::MIME"text/plain", model::AuroraModel)
    println(io, "AuroraModel:")
    println(io, "├── ", model.altitude_grid)
    println(io, "├── ", model.energy_grid)
    println(io, "├── ", model.pitch_angle_grid)
    println(io, "├── ", model.scattering)
    println(io, "├── ", model.ionosphere)
    println(io, "├── ", model.cross_sections)
    print(io, "└── B angle to zenith: $(model.B_angle_to_zenith)°")
end
