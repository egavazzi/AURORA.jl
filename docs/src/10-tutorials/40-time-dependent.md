# [Time-Dependent Simulation](@id Time-Dependent-Simulation)

This page expands on the flickering example from [Getting started](@ref) and focuses on
the time grid, `CFL_number`, and memory-aware loop partitioning.

## Setup

```@example time_dep
using AURORA

msis_file = find_msis_file()
iri_file  = find_iri_file()

model = AuroraModel(
    [100, 600],    # altitude limits [km]
    180:-30:0,     # pitch-angle bin edges [°]
    1000,          # maximum energy [eV]
    msis_file,
    iri_file,
    13             # magnetic field angle to zenith [°]
)
```

## Defining a flickering input flux

A sinusoidal flickering modulates the precipitating flux in time. The `z_source` parameter
introduces a realistic energy-dependent time delay — faster (higher energy) electrons
arrive before slower ones:

```@example time_dep
flux = InputFlux(
    FlatSpectrum(1e-2; E_min=100),        # 10 mW/m² above 100 eV
    SinusoidalFlickering(5.0);             # 5 Hz sinusoidal modulation
    beams=1,                               # single field-aligned beam
    z_source=3000.0                        # source at 3000 km altitude
)
```

Other modulation types include [`SquareFlickering`](@ref) (on/off square wave) and
[`SmoothOnset`](@ref) (smooth ramp-up). See the [Custom input flux](@ref "Custom Input Flux")
tutorial.

## Creating a time-dependent simulation

For time-dependent simulations, [`AuroraSimulation`](@ref) takes additional timing
parameters:

```@example time_dep
savedir = mkpath(joinpath("data", "time_dep_example"))
sim = AuroraSimulation(
    model,
    flux,
    0.2,                # total simulation time [s]
    0.01,               # output time step [s] (saves every 10 ms)
    savedir;            # save directory
    CFL_number=128,     # CFL stability parameter (higher → coarser internal stepping)
    max_memory_gb=4.0   # memory budget [GB] — controls loop partitioning
)
```

## Understanding the time grid

The `ResolvedTimeGrid` is constructed automatically inside `AuroraSimulation`. It
determines two key quantities:

- **`dt_resolved`**: the internal time step, chosen to satisfy the CFL condition. This is
  typically much smaller than the requested output `dt`.
- **`n_loop`**: the number of time chunks the simulation is split into. Each chunk is solved
  independently and its results are saved before moving to the next. This keeps memory
  usage within the specified `max_memory_gb` budget.

You can inspect these:

```@example time_dep
sim.time
```

If you only need the constant-input limit of this workflow, see
[Steady-State Simulation](@ref "Steady-State Simulation").

## Running the simulation

```@example time_dep
run!(sim)
nothing # hide
```

The solver loops over time chunks. Within each chunk, it sweeps through energies from high to
low (just like the steady-state case), but now advances the Crank-Nicolson scheme at each
energy to propagate the solution in time.

## Post-processing

Post-process the results as shown in [Getting started](@ref) Step 4.