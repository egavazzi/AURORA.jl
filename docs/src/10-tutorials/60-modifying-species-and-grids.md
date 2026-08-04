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

## Changing a species' density source

Each species carries a `density_source`: a callable mapping altitude (m) to density (m⁻³).
`initialize!` samples it onto the altitude grid. The built-in option is [`VectorDensity`](@ref)
(your own altitude/density vectors, pchip-interpolated in log-space). It is the universal way to
bring densities from any external atmospheric model (CCMC ModelWeb runs — MSIS 2.1, DTM,
WAM-IPE — radar inversions, …): reduce the source to two vectors and wrap it.

For a custom analytic profile, wrap it in the [`@law`](@ref) macro. `@law` captures the law's
source text so that it can be saved in `physics_state.jld2`. **Bare anonymous functions are rejected** 
as their source cannot be reconstructed when a saved model is reloaded.

```@example modular
# An analytic profile for N₂ (altitude in m → density in m⁻³)
model.species[:N2].density_source = @law h -> 1e18 .* exp.(-(h .- 100e3) ./ 30e3)
initialize!(model)
model.species[:N2].density
```

```julia
# Your own measured/downloaded profile, given as vectors (optionally tag its provenance):
model.species[:N2].density_source = VectorDensity(altitude_m, density_m3; source="ccmc_run.txt")

# Or read an existing AURORA-generated MSIS file (returns a VectorDensity):
model.species[:N2].density_source = MSISDensity(msis_file, :N2)
```

### Densities from another atmospheric model (e.g. a CCMC ModelWeb run)

For a whole atmosphere at once, [`NeutralProfile`](@ref) holds one [`VectorDensity`](@ref) per
species and can be passed straight to [`AuroraModel`](@ref) in place of an MSIS file. Two
readers build one for you:

```julia
neutrals = read_msis_file(msis_file)            # AURORA-generated MSIS file
neutrals = read_ccmc_msis("nrlmsis_output.txt") # CCMC ModelWeb NRLMSIS export

model = AuroraModel(altitude_lims, θ_lims, E_max, neutrals, electron)

# Index it to reach a single species
model.species[:N2].density_source = neutrals[:N2]
```

`read_ccmc_msis` handles the CCMC quirks: the variable-length preamble, densities in cm⁻³, and
the `9.999E-38` sentinel written where a species is not reported (those levels are dropped
per species).

For a model with no dedicated reader — DTM, WAM-IPE, a radar inversion — reduce the export to
two vectors and wrap them yourself:

```julia
using DelimitedFiles, AURORA

rows = readdlm("some_model_output.txt")
h_m  = Float64.(rows[2:end, 1]) .* 1e3     # km  → m
n_N2 = Float64.(rows[2:end, 2]) .* 1e6     # cm⁻³ → m⁻³

model.species[:N2].density_source = VectorDensity(h_m, n_N2; source="some_model_output.txt")
```

### When a law needs parameters, use a functor

`@law` is for **self-contained, closed-form** laws: the expression may reference only its
arguments and standard functions. A `@law` that closes over a local variable will be rejected
(because the captured value can't be saved as source). When a law must carry parameters or tabulated
data, define a small **functor** `struct` instead. Its fields are saved and can be reloaded with the
model. Here is an example:

```julia
struct ExponentialDensity      # parameters live in fields → can round-trip through JLD2
    n0::Float64                # density at z_ref (m⁻³)
    z_ref::Float64             # reference altitude (m)
    H::Float64                 # scale height (m)
end
(d::ExponentialDensity)(z) = d.n0 .* exp.(-(z .- d.z_ref) ./ d.H)

model.species[:N2].density_source = ExponentialDensity(1e18, 100e3, 30e3)
```

The same pattern applies to a parameterized cascading law or phase-function generator.
(If you dig into the source code you might find out that this is how the default 
atomic-oxygen cascading law is built in order to carry some tabulated coefficients)

## Changing the electron background (IRI)

The electron density and temperature (`ne`, `Te`) are the fifth argument to
[`AuroraModel`](@ref), supplied as an [`ElectronProfile`](@ref) — the electron analogue of
[`VectorDensity`](@ref). Like the density sources, it stores data (not a file path), so it
round-trips through `physics_state.jld2` and reproduces with no external file. Build one in any
of these ways:

```julia
# Run IRI-2020 for a date and position (via the bundled python iri2020; nothing written to disk):
electron = run_iri(; year=2005, month=10, day=8, hour=22, minute=0, lat=69.58, lon=19.23)

# From a CCMC ModelWeb IRI download:
electron = read_ccmc_iri("iri_output.txt")

# From an existing AURORA-generated IRI file (a bare path string also works as the 5th argument):
electron = read_iri_file(iri_file)

# Or your own vectors (altitude m, ne m⁻³, Te K):
electron = ElectronProfile(h_m, ne_m3, Te_K; source="my profile")

model = AuroraModel(altitude_lims, θ_lims, E_max, atmosphere, electron)
```

`read_ccmc_iri` handles the CCMC table directly (variable preamble, `Ne` in cm⁻³, `-1`
sentinels). For other electron-density sources, reduce them to `(h, ne, Te)` vectors and wrap
them in `ElectronProfile`, exactly as `VectorDensity` does for neutrals.

## Overriding cross-sections, phase functions, or cascading

The same interception window lets you replace other per-species physics before `initialize!`.
Set the *generator* (re-evaluated on every `initialize!`, so it tracks grid changes) rather than
the materialized array where possible:

```julia
# Phase-function generator: (θ, E) -> (elastic, inelastic) matrices.
# Can be a named function, functor, or @law, but not an anonymous function.
model.species[:N2].phase_fcn_generator = my_custom_phase_function

# Cascading: ionization thresholds + a secondary-electron distribution law f(E_s, E_p)
law = @law (E_s, E_p) -> 1/(15.2^2 + E_s^2)
model.species[:O2].cascading_spec = AURORA.CascadingSpec("O2", [12.07, 16.1], law)

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
law  = @law (E_s, E_p) -> 1.0 / (12.0^2 + E_s^2)  # we are completely inventing here
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
    Laws (density profiles, phase-function generators, cascading laws) must be reproducible so
    the model can be saved and reloaded. Use [`@law`](@ref) for closed-form laws, a functor
    `struct` when the law carries parameters, or a named function. **Bare anonymous functions
    are rejected**. The chosen law will be stored in `inputs/physics_state.jld2` and possible
    to restore on load.
