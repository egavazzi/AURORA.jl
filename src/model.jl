"""
    AuroraModel{AG, EG, PAG, IO, SP, FT, V}

Container for the grids, atmosphere, and collision data used by an AURORA simulation.

Construct cheaply with `AuroraModel(...)`, then call `initialize!(model)` (or `run!(sim)`)
to perform the heavy setup: scattering matrices, species densities, cross sections, and
cascading transfer matrices.
"""
mutable struct AuroraModel{AG<:AltitudeGrid, EG<:EnergyGrid, PAG<:PitchAngleGrid,
                            IO<:Ionosphere,
                            SP<:Tuple{Vararg{NeutralSpecies}},
                            FT, V<:AbstractVector{FT}}
    altitude_grid::AG
    energy_grid::EG
    pitch_angle_grid::PAG
    scattering::Union{Nothing, ScatteringData}   # nothing until initialize!(model)
    ionosphere::IO
    species::SP
    B_angle_to_zenith::FT
    s_field::V
    initialized::Bool
end

"""
    AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith=0)

Build a lightweight `AuroraModel`. Only the grids and ionosphere are set up; the expensive
computation (scattering matrices, species densities, cross sections, cascading data) is
deferred to `initialize!(model)`, which is called automatically by `run!(sim)`.

Between construction and initialization you can freely mutate species density profiles:
```julia
model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)
model.species[1].density_profile = VectorDensity(h, my_n2)
run!(sim)   # initialize!(model) picks up the custom profile here
```

## Inputs
- `altitude_lims`: altitude limits (km) for the bottom and top of the ionosphere
- `θ_lims`: pitch-angle beam limits (e.g. `180:-10:0`), 180° = field-aligned down
- `E_max`: upper energy-grid limit (eV)
- `msis_file`: path to the MSIS atmosphere file
- `iri_file`: path to the IRI ionosphere file
- `B_angle_to_zenith`: angle between the magnetic field and zenith (degrees, default 0)

## Returns
An uninitialized `AuroraModel`. Call `initialize!(model)` or `run!(sim)` to complete setup.
"""
function AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith=0)
    altitude_grid    = AltitudeGrid(altitude_lims[1], altitude_lims[2])
    energy_grid      = EnergyGrid(E_max)
    pitch_angle_grid = PitchAngleGrid(θ_lims)
    ionosphere       = Ionosphere(msis_file, iri_file, altitude_grid.h)
    species          = (N2Species(msis_file), O2Species(msis_file), OSpecies(msis_file))
    s_field          = altitude_grid.h ./ cosd(B_angle_to_zenith)

    FT = promote_type(eltype(s_field), typeof(B_angle_to_zenith))

    return AuroraModel{typeof(altitude_grid), typeof(energy_grid), typeof(pitch_angle_grid),
                       typeof(ionosphere), typeof(species),
                       FT, typeof(s_field)}(
        altitude_grid, energy_grid, pitch_angle_grid,
        nothing, ionosphere, species,
        FT(B_angle_to_zenith), s_field, false
    )
end

"""
    initialize!(model::AuroraModel; verbose=true, policy=CachePolicy())

Perform all heavy one-time setup for `model`:
1. Compute (or load from cache) the scattering matrices.
2. For each species: sample the density profile, load cross sections and excitation levels,
   build phase function matrices, and load/compute cascading transfer matrices.

Called automatically by `initialize!(sim)` → `run!(sim)`. Call it directly when you need
access to `model.scattering` or `model.species[i].density` before constructing a simulation.
"""
function initialize!(model::AuroraModel;
                     verbose::Bool = true,
                     policy::CachePolicy = CachePolicy())
    model.scattering = ScatteringData(model.pitch_angle_grid; verbose, policy)

    h  = model.altitude_grid.h
    eg = model.energy_grid
    θ  = model.scattering.θ_scatter

    for sp in model.species
        sp.density           = collect(Float64, sp.density_profile(h))
        apply_density_boundary!(sp.density)
        name_str             = String(sp.name)
        sp.cross_sections    = get_cross_section(name_str, eg.E_centers)
        sp.excitation_levels = load_excitation_threshold_for(name_str)
        sp.phase_fcn         = sp.phase_fcn_generator(θ, eg.E_centers)
        load_or_compute_cascading!(sp.cascading_data, eg; verbose, policy)
    end
    model.initialized = true
    return nothing
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
    println(io, "AuroraModel", model.initialized ? "" : " (not initialized)", ":")
    println(io, "├── ", model.altitude_grid)
    println(io, "├── ", model.energy_grid)
    println(io, "├── ", model.pitch_angle_grid)
    println(io, "├── ", isnothing(model.scattering) ? "ScatteringData: (not initialized)" : model.scattering)
    println(io, "├── ", model.ionosphere)
    println(io, "├── Species: ", join((String(sp.name) for sp in model.species), ", "))
    print(io, "└── B angle to zenith: $(model.B_angle_to_zenith)°")
end
