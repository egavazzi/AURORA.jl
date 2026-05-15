# [Architecture](@id Architecture)

This page gives a high-level overview of how AURORA is implemented: how the source code is
organized, the main types involved, and the execution flow from `run!` to the solver.

## Project layout

```
AURORA.jl/
├── data/                    # Simulation output (one subfolder per run)
│
├── docs/                    # Documentation (Documenter.jl)
│
├── ext/
│   ├── AURORA_viz.jl        # Makie extension (plotting & animation)
│   └── src/
│
├── internal_data/
│   ├── data_electron/       # IRI electron density/temperature files
│   ├── data_neutrals/       # MSIS neutral density/temperature files
│   ├── e_cascading/         # Cached cascading transfer matrices (per species)
│   │   ├── N2/
│   │   ├── O2/
│   │   └── O/
│   └── e_scattering/        # Cached scattering probability matrices
│
├── scripts/                 # Example scripts (plotting, analysis)
│
├── src/                     # Source code (see below)
│
└── test/                    # Unit tests
```

## Source code organization

```
src/
├── AURORA.jl                    # Module definition and exports
├── model.jl                     # AuroraModel constructor
│
├── grids/
│   ├── abstract_grid.jl         # AbstractGrid base type
│   ├── altitude_grid.jl         # AltitudeGrid
│   ├── energy_grid.jl           # EnergyGrid
│   └── pitch_angle_grid.jl      # PitchAngleGrid
│
├── ionosphere/
│   ├── ionosphere.jl            # Ionosphere struct (densities, temperatures)
│   ├── iri/                     # IRI model interface (Python iri2020)
│   └── msis/                    # MSIS model interface (Python pymsis)
│
├── input/
│   ├── spectra.jl               # Spectrum types (Flat, Gaussian, Maxwellian, File)
│   ├── modulations.jl           # Modulation types (Constant, Sinusoidal, Square, SmoothOnset)
│   └── input_flux.jl            # InputFlux, compute_flux
│
├── physics/
│   ├── cross_sections/          # e-N₂, e-O₂, e-O cross-section data
│   ├── cascading.jl             # Ionization cascading transfer matrices
│   ├── energy_degradation.jl    # Loss frequencies, scattering, source terms
│   ├── phase_functions.jl       # Differential cross sections → 3D scattering
│   └── scattering.jl            # Pitch-angle scattering matrices
│
├── solvers/
│   ├── types.jl                 # AbstractMode, SteadyStateMode, TimeDependentMode
│   ├── transport_matrices.jl    # TransportMatrices struct
│   ├── matrix_building.jl       # update_A!, update_B! (collision operators)
│   ├── sparse_indexing.jl       # Shared sparse infrastructure: BlockIndices, OperatorDiagonals, sparsity builders
│   ├── steady_state.jl          # Steady-state solver (in-place, allocation-free)
│   └── crank_nicolson.jl        # Crank-Nicolson time-dependent solver (in-place, allocation-free)
│
├── simulation/
│   ├── types.jl                 # AuroraSimulation, AbstractTimeConfig, RefinedTimeGrid, UniformTimeGrid
│   ├── cache.jl                 # SimulationCache, SolverCache, DegradationCache
│   ├── initialize.jl            # initialize! — allocates cache
│   └── run.jl                   # run!, solve!, energy_loop!, solve_energy_step!
│
├── analysis/
│   ├── analysis_types.jl        # Result types (VolumeExcitationResult, ColumnExcitationResult, IeTopResult)
│   ├── utilities.jl             # Helpers
│   ├── emissions.jl             # Volume and column excitation rates
│   ├── fluxes.jl                # Electron flux post-processing and downsampling
│   ├── heating.jl               # Electron heating rate
│   └── psd.jl                   # Phase space density analysis
│
├── utilities.jl                 # Helpers (v_of_E, CFL_criteria, beam_weight, ...)
└── precompiles.jl               # Precompilation workload
```

## Key types

| Type | Role |
|------|------|
| [`AuroraModel`](@ref) | Physical model: grids + atmosphere + cross sections |
| [`InputFlux`](@ref) | Precipitating electron specification: spectrum x modulation |
| [`AuroraSimulation`](@ref) | Complete simulation descriptor: model + flux + mode + output |
| [`AbstractMode`](@ref) | Solver strategy: [`SteadyStateMode`](@ref) or [`TimeDependentMode`](@ref) |
| `SimulationCache` | Internal workspace: solver matrices, flux arrays, factorizations |
| `SolverCache` | Per-energy sparse matrices `Mlhs`/`Mrhs`, KLU factorizations, `BlockIndices`, `OperatorDiagonals` |
| `TransportMatrices` | Matrices A (loss), B (scattering), D (diffusion), Q (source) |
| `BlockIndices` | Pre-computed `nzval` index arrays for a single block (replaces `Dict`-based mapping) |
| `OperatorDiagonals` | Dense diagonals of `Ddz_Up`, `Ddz_Down`, `Ddiffusion` (extracted once, reused each energy step) |
| [`RefinedTimeGrid`](@ref) | Time discretization with CFL refinement and loop partitioning. Last loop may have fewer save points 
than earlier loops. |
| [`UniformTimeGrid`](@ref) | Simple uniform grid for multi-step steady-state simulations |

## Execution flow

The entry point for all simulations is [`run!`](@ref). Here is the complete call chain:

```
run!(sim::AuroraSimulation)
│
├── initialize!(sim)                           # Allocate cache (if not done yet)
│   ├── Build SimulationCache
│   │   ├── Compute phase functions (once)
│   │   ├── Compute beam-to-beam scattering geometry
│   │   └── Load/compute cascading transfer matrices
│   └── Compute input flux → Ie_top array
│
├── Save parameters, neutral densities, Ie_top to disk
│
└── solve!(sim)
    │
    ├─── [steady-state: SteadyStateMode, single-step]
    │    └── energy_loop!(sim, Ie_top, 1, 1)
    │
    ├─── [steady-state: SteadyStateMode, multi-step]
    │    └── for i_step in 1:n_steps
    │        ├── Reset solver state
    │        ├── Compute Ie_top for this time point
    │        ├── energy_loop!(sim, Ie_top_step, 1, 1)
    │        └── Accumulate results
    │
    └─── [time-dependent: TimeDependentMode]
         └── for i_loop in 1:n_loop
             ├── Extract Ie_top chunk for this time window
             ├── energy_loop!(sim, Ie_top_chunk, i_loop, n_loop)
             └── Save results to disk


energy_loop!(sim, Ie_top, i_loop, n_loop)
│   Iterates over energies from HIGH to LOW (descending)
│
└── for iE in n_E:-1:1
    ├── update_matrices!(sim, iE)
    │   ├── update_A!()    # Loss frequencies from collisions
    │   ├── update_B!()    # Pitch-angle scattering probabilities
    │   ├── update_D!()    # Diffusion coefficients
    │   └── update_Ddiffusion!()  # Spatial diffusion operator
    │
    ├── solve_energy_step!(sim, iE, Ie_top)
    │   ├── [steady-state] → steady_state_scheme!()
    │   └── [time-dependent] → Crank_Nicolson!()
    │
    └── update_Q!(sim, iE)   # Source terms for lower energies
        ├── e-e collision energy transfer
        ├── Inelastic scattering contributions
        └── Ionization fragments (primary + secondary electrons)
```

**Note that the energy loop runs in descending order.** Because when an electron at energy
``E`` undergoes ionization, it produces secondary electrons at *lower* energies. By solving
from high to low, these secondary electron sources (stored in the Q matrix) are available
when the solver reaches the lower energy levels. This is what we refer to as the 
*cascading/degradation* in energy.

## Steady-state vs. time-dependent

The dispatch between steady-state and time-dependent is controlled by the
[`AbstractMode`](@ref) type stored in `sim.mode`:

- **Steady-state** ([`SteadyStateMode`](@ref)): No time stepping. The solver finds the
  equilibrium flux for each energy level by inverting a sparse linear system
  (``M \cdot I_e = Q``). With time parameters (`SteadyStateMode(duration=..., dt=...)`), each
  time point is solved independently (multi-step steady-state).

- **Time-dependent** ([`TimeDependentMode`](@ref)): The Crank-Nicolson scheme advances
  the solution in time at each energy level. The simulation may be split into multiple
  *loops* (time chunks) to fit within a memory budget, controlled by `max_memory_gb`.

For more details on the transport equation, matrix construction, and solver internals, see
[Solver Internals](@ref "Internals").
