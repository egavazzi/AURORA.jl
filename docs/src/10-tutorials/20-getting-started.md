# [Getting started](@id Getting-started)

This tutorial walks you through your first AURORA simulation—from setting up the
atmospheric model to running a short time-dependent simulation with a flickering input
flux.

## Loading AURORA

```@example getting_started
using AURORA
```

## Step 1: Create the atmospheric model

An [`AuroraModel`](@ref) bundles the numerical grids (altitude, energy, pitch angle),
the background ionosphere (neutral and electron densities from MSIS/IRI), and
collision cross-section data.

```@example getting_started
# Find (or generate) MSIS and IRI data files
# Default parameters: VISIONS-2 rocket launch conditions
# (2018-12-07 11:15 UTC, 76°N 5°E)
msis_file = find_msis_file()
iri_file  = find_iri_file()

model = AuroraModel(
    [100, 600],    # altitude limits [km]
    180:-45:0,     # pitch-angle bin edges [°] → 4 beams
    1000,          # maximum energy [eV]
    msis_file,
    iri_file,
    13             # magnetic field angle to zenith [°]
)
```

The constructor performs several steps automatically:
- Builds an [`AltitudeGrid`](@ref) with fine spacing (~150 m) at the lowest altitudes and 
  coarser spacing above.
- Builds an [`EnergyGrid`](@ref) with adaptive bin widths (finer at low energies 
  where numerous elastic and inelastic collisions need to be properly resolved).
- Builds a [`PitchAngleGrid`](@ref) from the given angle edges.
- Loads the neutrals (N₂, O₂, O) densities and temperature from the MSIS file.
- Loads electron densities and temperature from the IRI file.
- Precomputes electron collision cross sections and scattering data.

## Step 2: Define the precipitating electron flux

The precipitating electrons are specified by an [`InputFlux`](@ref), which combines an
*energy spectrum* with a *temporal modulation*.

For this first example we use a simple, 5 Hz flickering input.

```@example getting_started
flux = InputFlux(
  FlatSpectrum(1e-3; E_min=100),
  SinusoidalFlickering(5.0);
  beams=1,
  z_source=3000.0,
)
```

**Energy flux and modulation:**

- [`FlatSpectrum`](@ref)`(1e-3; E_min=100)` specifies a **total energy flux of 1 mW/m²**, uniform in energy above 100 eV. 
  The first argument is always the energy flux in W/m² (per square meter, per second).
- [`SinusoidalFlickering`](@ref)`(5.0)` applies a 5 Hz sinusoidal modulation with default
  amplitude 1.0, so the flux oscillates between **0 and 1 mW/m²**. If you want
  only 50% modulation depth (between 0.5 and 1 mW/m²), use `SinusoidalFlickering(5.0; amplitude=0.5)`.

**Beam and launch arguments:**

- The `beams=1` argument means that only the most field-aligned downward-going beam carries
  the precipitating flux.
- The `z_source` argument sets the altitude from which the electrons are launched. This
  introduces an energy-dependent travel-time delay, so higher-energy electrons reach the top
  of the ionosphere before lower-energy ones. 
- See the [Input flux](@ref "Input Flux")
  tutorial for other spectrum types, modulation types, and more details.

## Step 3: Create and run the simulation

Bundle the model and input flux into an [`AuroraSimulation`](@ref). For a time-dependent
simulation you specify the total duration and the output cadence:

```@example getting_started
savedir = mkpath(joinpath("data", "my_first_simulation"))
t = 0.2     # total simulation time [s]
dt = 0.01   # output time step [s] (save every 10 ms)
sim = AuroraSimulation(model, flux, t, dt, savedir; CFL_number=128)
```

Now run it:

```@example getting_started
run!(sim)
nothing # hide
```

Internally, AURORA refines the time grid as needed to satisfy the CFL criterion as given by
`CFL_number` (default `= 64`) and then sweeps through energies from high to low, advancing 
the Crank-Nicolson scheme for each energy.

## Step 4: Post-process

AURORA provides several analysis functions that compute derived quantities from the raw
electron flux and save them alongside the simulation output:

```@example getting_started
make_Ie_top_file(sim)              # boundary condition (input flux applied at top)
make_volume_excitation_file(sim)   # volumetric excitation rates for optical emissions
make_column_excitation_file(sim)   # column-integrated excitation rates
make_current_file(sim)             # field-aligned electron currents and energy fluxes
make_heating_rate_file(sim)        # electron heating rates
make_psd_file(sim)                 # electron phase-space density f(E, θ) and F(v∥)
```

The output files are saved in `data/my_first_simulation/`:

```@example getting_started
readdir(savedir)
```

## Next steps

- [Time-dependent simulation](@ref "Time-Dependent Simulation") — a deeper guide to the
  time grid, `CFL_number`, and memory partitioning.
- [Steady-state simulation](@ref "Steady-State Simulation") — the constant-input special case.
- [Custom input flux](@ref "Custom Input Flux") — explore all spectrum and modulation types.
