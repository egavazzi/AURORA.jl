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
    AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith=0;
                species = (N2Species(msis_file), O2Species(msis_file), OSpecies(msis_file)))

Build a lightweight `AuroraModel`. Only the grids and ionosphere are set up; the expensive
computation (scattering matrices, species densities, cross sections, cascading data) is
deferred to `initialize!(model)`, which is called automatically by `run!(sim)`.

## Inputs
- `altitude_lims`: altitude limits (km) for the bottom and top of the ionosphere
- `θ_lims`: pitch-angle beam limits (e.g. `180:-10:0`), 180° = field-aligned down
- `E_max`: upper energy-grid limit (eV)
- `msis_file`: path to the MSIS atmosphere file
- `iri_file`: path to the IRI ionosphere file
- `B_angle_to_zenith`: angle between the magnetic field and zenith (degrees, default 0)

## Keyword Arguments
- `species`: tuple (or any iterable) of [`NeutralSpecies`](@ref) instances that the model
  should simulate. Defaults to the standard atmosphere: N₂, O₂, O from `msis_file`.
  Pass a custom tuple to add, remove, or replace species:
  ```julia
  # Only two species:
  model = AuroraModel(...; species = (O2Species(msis_file), OSpecies(msis_file)))

  # Four species — custom 4th gas with pre-populated cross sections:
  custom_sp = NeutralSpecies(:MyGas, my_density_profile;
                             cascading_spec = my_spec, phase_fcn_generator = phase_fcn_N2)
  model = AuroraModel(...; species = (N2Species(msis_file), O2Species(msis_file),
                                      OSpecies(msis_file),  custom_sp))
  # Interception window: pre-populate custom cross sections before run!/initialize!
  model.species[end].cross_sections    = my_sigma_matrix   # [n_levels × n_E]
  model.species[end].excitation_levels = my_levels_matrix  # [n_levels × 2]
  run!(sim)
  ```

## Returns
An uninitialized `AuroraModel`. Call `initialize!(model)` or `run!(sim)` to complete setup.
"""
function AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith=0;
                     species = (N2Species(msis_file), O2Species(msis_file), OSpecies(msis_file)))
    altitude_grid    = AltitudeGrid(altitude_lims[1], altitude_lims[2])
    energy_grid      = EnergyGrid(E_max)
    pitch_angle_grid = PitchAngleGrid(θ_lims)
    ionosphere       = Ionosphere(msis_file, iri_file, altitude_grid.h)
    species_tuple    = Tuple(species)
    s_field          = altitude_grid.h ./ cosd(B_angle_to_zenith)

    FT = promote_type(eltype(s_field), typeof(B_angle_to_zenith))

    return AuroraModel{typeof(altitude_grid), typeof(energy_grid), typeof(pitch_angle_grid),
                       typeof(ionosphere), typeof(species_tuple),
                       FT, typeof(s_field)}(
        altitude_grid, energy_grid, pitch_angle_grid,
        nothing, ionosphere, species_tuple,
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
        if isempty(sp.cross_sections)
            sp.cross_sections    = get_cross_section(name_str, eg.E_centers)
        end
        if isempty(sp.excitation_levels)
            sp.excitation_levels = load_excitation_threshold_for(name_str)
        end
        sp.phase_fcn         = sp.phase_fcn_generator(θ, eg.E_centers)
        load_or_compute_cascading!(sp.cascading_data, eg; verbose, policy)
    end
    model.initialized = true
    return nothing
end


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
