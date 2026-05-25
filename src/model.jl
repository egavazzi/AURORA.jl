"""
    AuroraModel{AG, EG, PAG, SD, IO, SP, FT, V}

Container for the grids, atmosphere, and collision data used by an AURORA simulation.
"""
struct AuroraModel{AG<:AltitudeGrid, EG<:EnergyGrid, PAG<:PitchAngleGrid,
                   SD<:ScatteringData, IO<:Ionosphere,
                   SP<:Tuple{Vararg{NeutralSpecies}},
                   FT, V<:AbstractVector{FT}}
    altitude_grid::AG
    energy_grid::EG
    pitch_angle_grid::PAG
    scattering::SD
    ionosphere::IO
    species::SP
    B_angle_to_zenith::FT
    s_field::V
end

"""
    AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith=0;
                verbose=true, cache_policy=CachePolicy())

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
- `cache_policy=CachePolicy()`: controls reuse and persistence of the internal cascading and scattering cache

## Returns
An `AuroraModel` with fields:
- `altitude_grid::AltitudeGrid`
- `energy_grid::EnergyGrid`
- `pitch_angle_grid::PitchAngleGrid`
- `scattering::ScatteringData`
- `ionosphere::Ionosphere`
- `species::Tuple{NeutralSpecies, ...}`
- `B_angle_to_zenith::Real`
- `s_field::Vector`
"""
function AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith=0;
                     verbose=true, cache_policy::CachePolicy = CachePolicy())
    altitude_grid = AltitudeGrid(altitude_lims[1], altitude_lims[2])
    energy_grid = EnergyGrid(E_max)
    pitch_angle_grid = PitchAngleGrid(θ_lims)
    scattering = ScatteringData(pitch_angle_grid; verbose, policy=cache_policy)
    ionosphere = Ionosphere(msis_file, iri_file, altitude_grid.h)
    species = (
        N2Species(altitude_grid, energy_grid, scattering, msis_file),
        O2Species(altitude_grid, energy_grid, scattering, msis_file),
        OSpecies(altitude_grid, energy_grid, scattering, msis_file),
    )
    s_field = altitude_grid.h ./ cosd(B_angle_to_zenith)

    FT = promote_type(eltype(s_field), typeof(B_angle_to_zenith))

    return AuroraModel{typeof(altitude_grid), typeof(energy_grid), typeof(pitch_angle_grid),
                       typeof(scattering), typeof(ionosphere),
                       typeof(species),
                       FT, typeof(s_field)}(
        altitude_grid, energy_grid, pitch_angle_grid,
        scattering, ionosphere, species,
        FT(B_angle_to_zenith), s_field
    )
end

"""
    n_neutrals(model::AuroraModel)

Return neutral densities as a NamedTuple `(; nN2, nO2, nO)` sourced from the model's species.
"""
n_neutrals(model::AuroraModel) = (;
    nN2 = model.species[1].density,
    nO2 = model.species[2].density,
    nO  = model.species[3].density,
)

function Base.show(io::IO, model::AuroraModel)
    print(io, "AuroraModel($(model.altitude_grid), $(model.energy_grid), $(model.pitch_angle_grid))")
end

function Base.show(io::IO, ::MIME"text/plain", model::AuroraModel)
    println(io, "AuroraModel:")
    println(io, "├── ", model.altitude_grid)
    println(io, "├── ", model.energy_grid)
    println(io, "├── ", model.pitch_angle_grid)
    println(io, "├── ", model.scattering)
    println(io, "├── ", model.ionosphere)
    println(io, "├── Species: ", join((String(sp.name) for sp in model.species), ", "))
    print(io, "└── B angle to zenith: $(model.B_angle_to_zenith)°")
end
