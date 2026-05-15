# [Discretization](@id Discretization)

AURORA uses four discretizations that together determine what the model resolves and how
expensive it is to run:

- the [`AltitudeGrid`](@ref) which resolves vertical structure in the ionosphere,
- the [`EnergyGrid`](@ref) which resolves the electron energy distribution,
- the [`PitchAngleGrid`](@ref) which resolves the angular distribution into beams,
- and the simulation time grid which controls saved output times, and also 
    the finer internal solver steps in time-dependent runs.

The first three are part of the [`AuroraModel`](@ref) object. The time grid belongs to the
[`AuroraSimulation`](@ref) and depends on whether you run a single-step steady-state,
a multi-step steady-state, or a fully time-dependent simulation.

The resolution of the grids affect how well we capture different physics:

- altitude resolution controls how well sharp vertical structures are captured
- energy resolution controls how smoothly energy loss and cascading are represented
- pitch-angle resolution controls how well anisotropy and scattering are represented
- time resolution controls saved cadence, internal stability and how accuratly rapid changes are resolved

Increasing the resolution/size of these grids will improve the accuracy/fidelity of the 
simulations, but will also lead to longer runtimes. A balance between the two need to be found.


## Altitude grid

The altitude grid is built from the altitude limits passed to [`AuroraModel`](@ref), for
example `AuroraModel([100, 500], ...)`. The default internal grid is constructed with:

- uniform 150 m spacing below 100 km,
- nonuniform spacing above 100 km, with coarser and coarser steps as altitude increases
- but with a maximum spacing capped by `dz_max` (25 km by default).

The spacing can be inspected with `model.altitude_grid.Δh`.


## Energy grid

The energy grid is controlled by the `E_max` argument passed to [`AuroraModel`](@ref).
The default internal grid is piecewise non-uniform, with:

- a constant step size (~ 11.65 eV) above 500 eV to resolve ionization collisions
- finer and finer step sizes under 500 eV, down to 0.15 eV at 2eV, to resolve variations around
the inelastic collisions thresholds

The energy grid is defined by:

- `E_edges`: bin edges
- `E_centers`: bin centers
- `ΔE`: bin widths

and can inspected from `model.energy_grid`.

The requested `E_max` will be used as a constrain on the final bin center, i.e.
`E_centers[end]` will be always ≤ `E_max`. This means that the last bin edge `E_edges[end]` 
can lie above `E_max`.


## Pitch-angle grid

The pitch-angle grid construction is controlled by `θ_lims`, for example `180:-10:0`.

These values are beam edges, not beam centers. AURORA requires:

- `θ_lims` to include `180` and `0`
- and the values to be provided in descending order

The transport equation is written in terms of the pitch-angle cosine `μ = cos(θ)`, so the
grid also stores:

- `μ_lims`, the cosine of the angle limits
- `μ_center`, the effective beam centers used by the transport solver

These are available from `model.pitch_angle_grid`.


## Time grid

The time grid is part of the simulation setup rather than the model discretization itself.

There are three relevant cases:

- [`SingleStepConfig`](@ref): one steady-state solve, no time evolution
- [`UniformTimeGrid`](@ref): multiple saved times in steady-state mode, where each time point
  is solved independently
- [`RefinedTimeGrid`](@ref): time-dependent mode, where the save cadence is internally
  refined to satisfy the CFL condition

For time-dependent simulations, `sim.time.dt` is the saved output cadence. The solver
then derives a smaller `dt_internal` based on the CFL condition and solves over it. 
The simulation will be partitioned into several loops (`n_loop`) if needed so the
working arrays stay within the memory (RAM) budget.

Expect `dt_internal` to become much smaller than `dt` in time-dependent runs, especially
for runs with a high `E_max` energy (faster electrons motion will require a finer internal time grid).

The time grid configuration can be inspected from the simulation object `sim.time`.


## Matching external inputs to the model grids

When loading external fluxes from a file, the file dimensions must match the model grids used
by the simulation. In particular, AURORA does not interpolate in pitch angle or energy when
reading file-based fluxes, so the file must be prepared on the same beam and energy grid as
the target model.

The safest workflow is:

1. build the [`AuroraModel`](@ref) with the settings you plan to simulate,
2. read `model.pitch_angle_grid.μ_center` and `model.energy_grid.E_centers`,
3. write the external file on those exact grids.

See the [Input Flux](@ref "Input Flux") tutorial for the concrete file format and loading
workflow.