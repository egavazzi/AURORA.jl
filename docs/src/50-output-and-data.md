# [Output & data](@id Output)

```@meta
CurrentModule = AURORA
```

[`run!(sim)`](@ref run!) writes all of a simulation's output into its save directory
(`sim.output.savedir`). Nothing is written during construction or [`initialize!`](@ref). 
The format is split by job:

| Format | Used for |
|--------|----------|
| **NetCDF** (`.nc`) | N-dimensional array data (fluxes, atmosphere, analysis)    |
| **JLD2** (`.jld2`) | the full Julia `AuroraModel` (species + physics matrices)  |
| **TOML** (`.toml`) | scalar parameters, for quick inspection                    | 

## Directory layout

```
savedir/
├── config.toml               # scalar parameters (grid limits, dt, E_max, version, commit)
├── inputs/
│   ├── atmosphere.nc         # neutral & electron density / temperature profiles
│   └── physics_state.jld2    # full AuroraModel (species + physics matrices) — Julia reload
├── simulation_data.nc        # electron flux output (appended per solver loop, crash-safe)
└── analysis/                 # derived quantities, written by the make_* functions
    ├── volume_excitation.nc
    ├── column_excitation.nc
    ├── Ie_top.nc
    ├── currents.nc
    ├── heating_rate.nc
    └── psd.nc
```

`config.toml`, `inputs/`, and `simulation_data.nc` are written by `run!`. The `analysis/`
files are written separately by the post-processing functions (see [Analysis](@ref)).

## `simulation_data.nc`

The main output. The `time` dimension is appended to after each solver loop.

| Variable | Dims | Units | Notes |
|----------|------|-------|-------|
| `altitude` | `(altitude,)` | m | altitude grid centers |
| `pitch_angle` | `(pitch_angle,)` | 1 | cosine of pitch angle, beam center (μ) |
| `energy` | `(energy,)` | eV | energy bin centers |
| `time` | `(time,)` | s | save-cadence time axis |
| `energy_edges` | `(energy_bounds,)` | eV | energy bin edges |
| `mu_lims` | `(pitch_angle_bounds,)` | 1 | μ bin boundaries |
| `dE` | `(energy,)` | eV | energy bin width |
| `beam_weight` | `(pitch_angle,)` | sr | solid-angle beam weight |
| `Ie` | `(altitude, pitch_angle, time, energy)` | m⁻² s⁻¹ | electron **number** flux |
| `Ie_input` | `(pitch_angle, time_input, energy)` | m⁻² s⁻¹ | input precipitation (boundary condition) |
| `time_input` | `(time_input,)` | s | time axis for `Ie_input` |

Global attributes: `aurora_version`, `commit_hash`, `creation_time`.

!!! warning "`Ie` is a number flux, not a differential flux"
    `Ie` is already integrated over each energy bin **and** over each beam's solid angle, so
    its units are m⁻² s⁻¹ (electrons per m² per second, per bin, per beam) — *not*
    eV⁻¹ m⁻² s⁻¹ sr⁻¹. To obtain a differential flux, divide by `dE` and
    `beam_weight` (both saved in the file).

!!! note "`Ie_input` vs `analysis/Ie_top.nc`"
    `Ie_input` (here) is the **boundary condition** fed into the simulation. `analysis/Ie_top.nc`
    (written by [`make_Ie_top_file`](@ref)) is the simulation **output** at the top altitude
    slice — it contains *all* beams, including the upward (backscattered) flux. They are
    different quantities; see [`load_Ie_top`](@ref).

## `inputs/`

- `atmosphere.nc` — altitude grid, electron density `ne` and temperature `Te`, and one number
  density variable per neutral species (`nN2`, `nO2`, `nO`, …).
- `physics_state.jld2` — the complete [`AuroraModel`](@ref), including the materialized
  scattering and cascading matrices. Reload it in Julia with
  `model = JLD2.load("savedir/inputs/physics_state.jld2", "model")`.

## Controlling output — [`AuroraOutputManager`](@ref)

`AuroraSimulation` takes an [`AuroraOutputManager`](@ref) as its third argument, which holds
all output options:

```julia
out = AuroraOutputManager("my_run"; overwrite = false, compress = true)
sim = AuroraSimulation(model, flux, out; mode = TimeDependentMode(duration = 0.5, dt = 0.001))

# Convenience: passing a plain String wraps it in an AuroraOutputManager with the defaults
sim = AuroraSimulation(model, flux, "my_run"; mode = SteadyStateMode())
```

- `overwrite = false` — `run!` errors if `simulation_data.nc` already exists in `savedir`
  (pass `overwrite = true` to allow replacing it).
- `compress = true` — zlib compression (`deflatelevel = 4`) on all NetCDF variables.
- `savedir` may be a relative or absolute path and is created by `run!`. An empty or
  whitespace `savedir` falls back to `backup/<yyyymmdd-HHMM>/` in the current working directory.

## Reading results

In Julia, use the loader helpers (see [Analysis](@ref) for the result types):

```julia
result = load_results("my_run")                # SimulationResult: Ie, t, h_atm, E_centers, E_edges, dE, mu_lims
vol    = load_volume_excitation("my_run")      # after make_volume_excitation_file
top    = load_Ie_top("my_run")                 # after make_Ie_top_file
```

Because the files are self-describing, no AURORA-specific code is needed elsewhere. In Python:

```python
import xarray as xr
ds = xr.open_dataset("my_run/simulation_data.nc")
Ie = ds["Ie"]                                  # dims: altitude, pitch_angle, time, energy
profile = Ie.sel(time=0.1, method="nearest").sum("pitch_angle")
```

In MATLAB:

```matlab
Ie  = ncread("my_run/simulation_data.nc", "Ie");
alt = ncread("my_run/simulation_data.nc", "altitude");
```

## Analysis outputs

The `analysis/` files are produced on demand by the post-processing functions, each reading
`simulation_data.nc` (and, where needed, `inputs/atmosphere.nc`). See
[Post-Processing & Analysis](@ref Post-Processing) for usage and the [Analysis](@ref) API page
for the per-function compatibility table.
