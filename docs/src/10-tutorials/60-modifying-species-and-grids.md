# Modifying species and grids

An [`AuroraModel`](@ref) is built in **two phases**:

1. `AuroraModel(...)` is *lightweight* — it sets up the grids and reads the background
   atmosphere, but does **not** compute the expensive data (densities sampled on the grid,
   cross-sections, phase functions, scattering and cascading matrices).
2. [`initialize!`](@ref) does that heavy work. It is called automatically by `run!(sim)`.

The gap between the two is an **interception window**: anything you change on the model before
`initialize!`/`run!` is picked up when the derived data is built. This is how you can customize
the atmosphere without needing to put your hands in the internals.

## Changing the grids

The grid fields can be reassigned on an existing model. Doing so automatically marks the model
as uninitialized, so the next `run!(sim)` rebuilds everything for the new grid — no manual
re-initialization needed.

```@example modular
using AURORA
msis_file = find_msis_file()
iri_file  = find_iri_file()

model = AuroraModel([100, 300], 180:-90:0, 100, msis_file, iri_file)
initialize!(model)
model.initialized
```

```@example modular
model.energy_grid = EnergyGrid(500)   # also: model.altitude_grid, model.pitch_angle_grid
model.initialized                     # → false: derived data is now stale
```

```@example modular
initialize!(model)                    # rebuilds densities, cross-sections, … for the new grid
model.energy_grid.n
```

When a simulation already exists, just change the grid and call `run!(sim)` — it detects the
change and rebuilds the model and its cache before solving.

## Changing a species' density profile

Each species carries a `density_profile`: any callable mapping altitude (m) to density (m⁻³).
`initialize!` samples it onto the altitude grid. Built-in options are [`MSISDensity`](@ref) and
[`VectorDensity`](@ref) (your own altitude/density vectors, that will get pchip interpolated
in log-space), but any function works.

```@example modular
# A plain analytic profile for N₂ (altitude in m → density in m⁻³)
model.species[:N2].density_profile = h -> 1e18 .* exp.(-(h .- 100e3) ./ 30e3)
initialize!(model)
model.species[:N2].density[:N2]
```

```julia
# Your own measured/downloaded profile, given as vectors:
model.species[:N2].density_profile = VectorDensity(altitude_m, density_m3)

# Or the default MSIS-backed profile, found/created by find_msis_file():
model.species[:N2].density_profile = MSISDensity(msis_file, :N2)
```

## Overriding cross-sections, phase functions, or cascading

The same interception window lets you replace other per-species physics before `initialize!`.
Set the *generator* (re-evaluated on every `initialize!`, so it tracks grid changes) rather than
the materialized array where possible:

```julia
# Phase-function generator: (θ, E) -> (elastic, inelastic) matrices
model.species[:N2].phase_fcn_generator = my_phase_function

# Cascading: ionization thresholds + a secondary-electron distribution law f(E_s, E_p)
model.species[:O2].cascading_spec = AURORA.CascadingSpec("O2", [12.07, 16.1], (E_s, E_p) -> 1/(15.2^2 + E_s^2))

# Cross-sections can be pre-populated directly for a non-standard species (see below)
```

## Adding, removing, or replacing species

Pass an explicit `species` tuple to the constructor. The defaults are
`N2Species`/`O2Species`/`OSpecies`, which are helper functions accepting
a density source (i.e. msis_file).

```julia
# Two species only:
model = AuroraModel(alt_lims, θ_lims, E_max, msis_file, iri_file;
                    species = (O2Species(msis_file), OSpecies(msis_file)))
```

A completely custom species needs its cascading law and a phase-function generator. Because the
built-in cross-section library only knows N₂/O₂/O, pre-populate the cross-sections and
excitation levels for a new gas in the interception window:

```julia
law  = (E_s, E_p) -> 1.0 / (12.0^2 + E_s^2)  # we are completely inventing here
spec = AURORA.CascadingSpec("Ar", [15.76, 27.63], law)
argon = NeutralSpecies(:Ar, MSISDensity(msis_file, :Ar);
                       cascading_spec = spec, phase_fcn_generator = phase_fcn_N2)

model = AuroraModel(alt_lims, θ_lims, E_max, msis_file, iri_file;
                    species = (N2Species(msis_file), O2Species(msis_file),
                               OSpecies(msis_file), argon))

model.species[:Ar].cross_sections    = my_sigma_matrix   # [n_levels × n_E]
model.species[:Ar].excitation_levels = my_levels_matrix  # [n_levels × 2]

run!(AuroraSimulation(model, flux, savedir; mode))
```

!!! tip
    Define custom profiles, phase functions, and cascading laws as small `struct`s (functors) or
    named functions rather than anonymous closures. In-line functions do work, but only 
    proper structs/named functions can be cleanly saved to file. To be fair we do not save
    them to file right now, but this is in the plans for improved reproducibility.
