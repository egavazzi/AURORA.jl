# [Time-Dependent Simulation](@id Time-Dependent-Simulation)

This tutorial walks you through your first AURORA simulation—from setting up the model to running a short time-dependent simulation with a flickering input flux.

## Loading AURORA

```@example time_dep
using AURORA
```

## Step 1: Create the atmospheric model

An [`AuroraModel`](@ref) bundles the numerical grids (altitude, energy, pitch angle),
the background ionosphere (neutral and electron densities from MSIS/IRI), and
collision cross-section data.

```@example time_dep
# Find (or generate) MSIS and IRI data files
# Default parameters: VISIONS-2 rocket launch conditions
# (2018-12-07 11:15 UTC, 76°N 5°E)
msis_file = find_msis_file()
iri_file  = find_iri_file()

model = AuroraModel(
    [100, 600],    # altitude limits [km]
    180:-30:0,     # pitch-angle bin edges [°] → 6 beams
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

For this first example we use a simple 5 Hz flickering input.

```@example time_dep
flux = InputFlux(
  FlatSpectrum(1e-3; E_min=100),
  SinusoidalFlickering(5.0);
  beams=1,
  z_source=3000.0,
)
```

**Energy flux and modulation:**

- [`FlatSpectrum`](@ref)`(1e-3; E_min=100)` specifies a **total energy flux of 1 mW/m²**, uniform in energy above 100 eV. 
- [`SinusoidalFlickering`](@ref)`(5.0)` applies a 5 Hz sinusoidal modulation with default
  amplitude 1.0, so the flux oscillates between **0 and 1 mW/m²**.

See the [Input flux](@ref "Input Flux") tutorial for other spectrum types, modulation types, and more details.

**Beam and launch arguments:**

- The `beams=1` argument means that only the most field-aligned downward-going beam carries
  the precipitating flux.
- The `z_source` argument sets the altitude from which the electrons are launched. This
  introduces an energy-dependent travel-time delay, so higher-energy electrons reach the top
  of the ionosphere before lower-energy ones. 

## Step 3: Create and run the simulation

Bundle the model and input flux into an [`AuroraSimulation`](@ref). For a time-dependent
simulation you specify the total duration and the output cadence:

```@example time_dep
savedir = mkpath(joinpath("data", "my_first_simulation"))
savedir = mktempdir()  # hide — redirect to OS temp so .mat files are not deployed to gh-pages

sim = AuroraSimulation(
    model,
    flux,
    0.2,                # total simulation time [s]
    0.01,               # output time step [s] (save every 10 ms)
    savedir;
    CFL_number=128,     # CFL stability parameter (higher → coarser internal stepping)
    max_memory_gb=4.0   # memory budget [GB] — controls loop partitioning
)
```

### Understanding the time grid

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

### Running

```@example time_dep
run!(sim)
nothing # hide
```

Internally, AURORA refines the time grid as needed to satisfy the CFL criterion and then
sweeps through energies from high to low, advancing the Crank-Nicolson scheme for each
energy. The solver loops over time chunks; within each chunk it solves the full energy
sweep and saves the results before moving to the next.

## Step 4: Post-process

AURORA provides several analysis functions that compute derived quantities from the raw
electron flux and save them alongside the simulation output:

```@example time_dep
make_Ie_top_file(sim)              # boundary condition (input flux applied at top)
make_volume_excitation_file(sim)   # volumetric excitation rates for optical emissions
make_column_excitation_file(sim)   # column-integrated excitation rates
make_current_file(sim)             # field-aligned electron currents and energy fluxes
make_heating_rate_file(sim)        # electron heating rates
make_psd_file(sim)                 # electron phase-space density f(E, θ) and F(v∥)
```

The output files are saved in `data/my_first_simulation/`:

```@example time_dep
readdir(savedir)
```

## Step 5: Visualize

AURORA.jl provides helper plotting functions to visualize results through a [Makie](https://github.com/MakieOrg/Makie.jl) extension. 
Install and load a Makie backend to access them (more information about this in [`Visualization`](@ref "Visualization")).

### Volume excitation rate

```@example time_dep
using CairoMakie
CairoMakie.activate!() # hide
vol = load_volume_excitation(savedir)

fig = Figure()
ax  = Axis(fig[1, 1]; xlabel = "Time (s)", ylabel = "Altitude (km)", title = "4278 Å")
hm  = plot_excitation!(ax, vol; field = :Q4278)
Colorbar(fig[1, 2], hm; label = "photons/m³/s")
fig
```

### Column emission intensity

```@example time_dep
using CairoMakie
CairoMakie.activate!() # hide
col = load_column_excitation(savedir)

fig = Figure()
ax  = Axis(fig[1, 1]; xlabel = "Time (s)", ylabel = "Intensity (R)", yscale = log10)
plot_column_excitation!(ax, col)
Legend(fig[1, 2], ax)
fig
```

### Animation

[`animate_Ie_in_time`](@ref) produces an animation of the electron flux distribution
plotted over altitude and energy, with different panels for different pitch-angle streams,
stepping over time:

```@example time_dep
using CairoMakie
CairoMakie.activate!() # hide
animate_Ie_in_time(savedir; framerate = 15)
nothing # hide
```

```@raw html
<video autoplay loop muted playsinline controls src="../data/my_first_simulation/animation.mp4"/>
```

For the full API of all visualization functions, see [Visualization](@ref).

## Next steps

- [Steady-state simulation](@ref "Steady-State Simulation") — the steady-state special case.
- [Input flux](@ref "Input Flux") — explore all spectrum and modulation types.
- [Post-processing & analysis](@ref "Post-Processing") — detailed walkthrough of all analysis functions.