# Folder structure

The package is structured as follows:

```
AURORA.jl/
├── data/
│   └── <simulation_name>/...
│   └── ...
│
├── docs/...
│
├── internal_data/
│   ├── data_electron/
│   ├── data_neutrals/
│   ├── e_cascading/
│   │   ├── N2/
│   │   ├── O/
│   │   └── O2/
│   └── e_scattering/
│
├── scripts/
│   ├── Control_script.jl
│   ├── Control_script_SteadyState.jl
│   ├── analysing/
│   │   ├── Calculate_densities.jl
│   │   └── Downsampling_fluxes.jl
│   └── plotting/
│       ├── Compare_Ie.jl
│       ├── Compare_emission_rate.jl
│       ├── Compare_ionization_rate.jl
│       ├── Control_animations.jl
│       ├── Plot_densities.jl
│       ├── Plot_emission_rate.jl
│       ├── Plot_ionization_rate.jl
│       ├── Plot_phase_functions.jl
│       └── Plot_setup.jl
│
├── src/
│   ├── AURORA.jl
│   ├── main.jl
│   ├── setup.jl
│   ├── input.jl
│   ├── cascading.jl
│   ├── scattering.jl
│   ├── phase_functions.jl
│   ├── energy_degradation.jl
│   ├── matrix_building.jl
│   ├── crank_nicolson.jl
│   ├── steady_state.jl
│   ├── analysis.jl
│   ├── utilities.jl
│   ├── iri/
│   └── msis/
│
└── test/
    ├── runtests.jl
    └── ...
```

The folder `data/` contains the subfolders where simulation results are saved. Each simulation run creates its own subfolder named after the simulation.

The folder `docs/` contains all the necessary scripts to power this documentation.

The folder `internal_data/` holds data used internally by the package:
- `data_electron/` and `data_neutrals/` contain cross-section and atmospheric species data.
- `e_cascading/` is where the electron cascading data produced by the simulations are cached for future use. The data are organized by species in the subfolders `N2/`, `O2/` and `O/`. You should not need to venture into this folder.
- `e_scattering/` contains precomputed electron scattering matrices.

The folder `scripts/` contains scripts for users to run simulations and analyse or plot results:
- `Control_script.jl` and `Control_script_SteadyState.jl` are the main entry points for running time-dependent and steady-state simulations respectively.
- `analysing/` contains scripts for post-processing simulation outputs.
- `plotting/` contains scripts for visualising simulation results.

The folder `src/` contains the source code of the model. The `iri/` and `msis/` subfolders contain the interfaces to the IRI and MSIS atmospheric models respectively.

The folder `test/` contains the automated tests for the package.
