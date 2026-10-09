# Changelog

## Unreleased
- **Breaking** A neutral species describes its collisions with an in-code channel table instead of the `internal_data/data_neutrals/<species>_levels.dat` / `.name` files [#186](https://github.com/egavazzi/AURORA.jl/pull/186)
  - New type `CollisionChannel(name, cross_section, energy_loss, n_secondaries; source)` for one inelastic channel; `CollisionChannel(c; energy_loss = ...)` copies it with fields replaced.
  - `NeutralSpecies` gains `elastic_cross_section`, `channels` and `secondary_law`; `cross_sections`, `excitation_levels`, `cascading_spec` and `cascading_data` are rebuilt from them at every `initialize!(model)` (row 1 elastic, row `i + 1` = `channels[i]`). Assigning one of these four fields throws; `run!` throws when the channel table was edited after `initialize!(model)`.
  - **Breaking** `NeutralSpecies(name, density_source; ...)` takes `elastic_cross_section`, `channels` and `secondary_law` in place of `cascading_spec`. `NeutralSpecies` is no longer parametric.
  - New functions `default_channels(:N2)`, `default_elastic_cross_section`, `default_secondary_law`, `channel_names` and `ionizing_channels`, and the unexported `AURORA.default_cascading_spec` and `AURORA.channel`. The built-in tables live in `src/physics/cross_sections/channels_N2.jl`, `channels_O2.jl` and `channels_O.jl`.
  - `AURORA.cascading_spec_from_channels(name, secondary_law, channels)` derives the ionization thresholds and secondary counts from a channel table, one per distinct (energy loss, secondary count) pair.
  - **Breaking** `DefaultCascadingSpecN2` / `O2` / `O`, `load_excitation_threshold`, `load_excitation_threshold_for` and `get_level_names` are removed. `get_cross_section` also accepts a `Symbol` species name.
  - `make_volume_excitation_file` takes the ionization cross sections from the run's `inputs/physics_state.jld2`, which it now requires.
  - The cross sections and excitation levels of the built-in N₂, O₂ and O species are unchanged, bit for bit.
- A cross section is now set to zero in the energy bins whose center is below the channel's energy loss, with a warning [#185](https://github.com/egavazzi/AURORA.jl/pull/185)
- **Numerical Breaking (small)** Non-ionizing degradation no longer renormalizes the partition over lower bins; electrons degraded below the lowest grid edge thermalise (are removed) instead of being pushed into the first bin [#185](https://github.com/egavazzi/AURORA.jl/pull/185)
- **Numerical Breaking (small)** Ionizing collisions near a threshold no longer delete the electron, and both cascading spectra are normalized by the row's ionization event count, so outgoing electrons below the lowest energy bin thermalise (are removed) instead of being redistributed on-grid [#185](https://github.com/egavazzi/AURORA.jl/pull/185)
- **Breaking** `initialize!(model)` throws an `ArgumentError` instead of warning when an energy bin is wider than a species' lowest ionization threshold [#184](https://github.com/egavazzi/AURORA.jl/pull/184)
- Fix the `e_N2*`, `e_O2*` and `e_O*` cross-section functions for unsorted or integer input energies [#184](https://github.com/egavazzi/AURORA.jl/pull/184)
- A cascading cache file now checks if its thresholds, secondary counts and secondary law match [#184](https://github.com/egavazzi/AURORA.jl/pull/184)
- Add an energy-budget diagnostic (analysis function) which reports how much of the precipitating energy flux goes into neutral excitation and ionization (split per channel and per species), thermal-electron heating, backscatter out of the top and absorption at the bottom of the grid, plus the unaccounted residual [#155](https://github.com/egavazzi/AURORA.jl/pull/155)
- **Breaking** `AuroraModel` takes the neutral atmosphere and electron background as data instead of MSIS/IRI file paths [#166](https://github.com/egavazzi/AURORA.jl/pull/166)
  - New types `NeutralAtmosphere` (one `DensityProfile` per species, indexed as `neutrals[:N2]`), `DensityProfile` and `ElectronProfile`. They hold the data itself, so a model saved to `physics_state.jld2` reloads without the original files, and carry a free-form `origin` string written into `inputs/atmosphere.nc`. File paths are still accepted, and read at construction.
  - New functions `run_msis`, `read_msis_file`, `read_ccmc_msis` (returning a `NeutralAtmosphere`) and `run_iri`, `read_iri_file`, `read_ccmc_iri` (returning an `ElectronProfile`).
  - `run_msis` and `run_iri` take a `save_to` directory in which to write the model output as an AURORA text file. `find_msis_file` / `find_iri_file` remain the cached route: they reuse a matching file from the package's file store, and compute and write one when there is none.
  - **Breaking** `NeutralSpecies.density_profile` is renamed `density_source`. `MSISDensity` and `VectorDensity` are removed; use `read_msis_file(file)[:N2]` or `DensityProfile(h, n)`.
  - **Breaking** `Ionosphere` is now `Ionosphere(electron_source, h_atm)` and no longer stores `msis_file`/`iri_file`.
  - Sampling a profile outside its native altitude range now warns that the values there are extrapolated.
  - Species that MSIS does not report at low altitude (N, anomalous O) no longer produce `NaN` densities: each species keeps only the levels where it is defined.
- **Breaking** Rename `AuroraSimulation.cache` to `AuroraSimulation.workspace`, and replace `cache_initialized` with `workspace.initialized`.
  The simulation working-state types are also renamed from cache to workspace, e.g. `SolverCache` → `SolverWorkspace`, `DegradationCache` → `DegradationWorkspace` [#161](https://github.com/egavazzi/AURORA.jl/pull/161)
- Cascading matrices built from an `@law` secondary law (e.g. default N₂, O₂) are now much faster to compute (~36x for a full 30 keV N₂ build, with ~200x less memory allocated), by avoiding dynamic dispatch [#175](https://github.com/egavazzi/AURORA.jl/pull/175)
- Compute the double-ionization cascading matrices with a numerical-CDF method with fixed Gauss–Legendre rules instead of adaptive 3-D cubature. Make it possible to use very large energy grids (> 100 keV) [#174](https://github.com/egavazzi/AURORA.jl/pull/174)
- Faster single-ionization cascading matrix calculations and better report progress [#169](https://github.com/egavazzi/AURORA.jl/pull/169)
- Remove the ad-hoc spatial diffusion operator (`D·∂²Ie/∂z²`) from both the steady-state and time-dependent solvers [#168](https://github.com/egavazzi/AURORA.jl/pull/168)
  - The operator was meant to model the spread in arrival times of electrons within a finite (E, μ) bin, but in its current form it contributed nothing measurable, and it did not belong in the steady-state equations in the first place. The numerical diffusion of the advection scheme already produces a comparable spread at default resolution.

## v0.8.0 - 2026-07-28
- **Breaking** :sparkles: New simulation interface :sparkles: [#114](https://github.com/egavazzi/AURORA.jl/pull/114) [#125](https://github.com/egavazzi/AURORA.jl/pull/125) [#126](https://github.com/egavazzi/AURORA.jl/pull/126) [#138](https://github.com/egavazzi/AURORA.jl/pull/138)
  - Simulations are now set up by building an `AuroraModel`, constructing an `AuroraSimulation`, and calling `run!(sim)`. The functions `calculate_e_transport()` and `calculate_e_transport_steady_state()` are removed.
  - Neutral species can be inspected, modified, removed or even added.
  - The grids, cross-section data and simulation state are now proper types, which makes them easier to inspect.
  - Visit the updated online [documentation](https://egavazzi.github.io/AURORA.jl/dev/) for more details and examples.
- **Breaking** New output format [#140](https://github.com/egavazzi/AURORA.jl/pull/140)
  - Results are saved as NetCDF/TOML/JLD2 instead of `.mat` files.
  - `savedir` is now an absolute path or a path relative to the current directory, instead of always being placed under the package `data/` folder.
  - The simulation model state is always saved to disk next to the results, and can be reloaded for full reproducibility.
- **Breaking** `IeE_tot` is now the field-aligned (vertical) energy flux entering the top of the ionosphere, instead of the omnidirectional energy flux [#150](https://github.com/egavazzi/AURORA.jl/pull/150)
  - Previously the normalization did not project the flux onto the field line, so the energy actually deposited in the atmosphere depended on the pitch-angle beams selected as input (as little as half of the requested `IeE_tot` for wide selections). `IeE_tot` now denotes the vertical energy flux that actually enters the atmosphere, and is invariant with respect to the beam selection.
  - Simulations from previous versions can be reproduced by rescaling `IeE_tot`, as the transport is linear in the magnitude of the input flux.
- **Numerical Breaking (small)** Minor correction of the Crank-Nicolson top boundary indexing/timing [#120](https://github.com/egavazzi/AURORA.jl/pull/120)
- **Numerical Breaking (small)** Refactor of the cascading functions, during which a missing factor was found and fixed [#130](https://github.com/egavazzi/AURORA.jl/pull/130)
- **Numerical Breaking (small)** The last `E_centers` is now always ≤ than `E_max` [#136](https://github.com/egavazzi/AURORA.jl/pull/136)
- **Numerical Breaking (small)** Remove the erf taper of density profiles at the top of the ionosphere [#141](https://github.com/egavazzi/AURORA.jl/pull/141)
- **Numerical Breaking (small)** Improve physical accuracy of cascading calculations [#132](https://github.com/egavazzi/AURORA.jl/pull/132)
- **Numerical Breaking (small)** Fix a bug in the energy cascading that led to a small creation of energy [#153](https://github.com/egavazzi/AURORA.jl/pull/153)
- **Numerical Breaking (small)** Properly handle double-ionization collisions [#156](https://github.com/egavazzi/AURORA.jl/pull/156)
- Add error message for when iri calculations return invalid data or when loading invalid data from file [#116](https://github.com/egavazzi/AURORA.jl/pull/116)
- Switch from IRI2016 to IRI2020 model, solving the issue with recent dates that could not be computed [#117](https://github.com/egavazzi/AURORA.jl/pull/117)
- Handle invalid iri values at top and bottom ends [#118](https://github.com/egavazzi/AURORA.jl/pull/118)
- Cached cascading and scattering matrices created with different versions of AURORA are now automatically skipped [#135](https://github.com/egavazzi/AURORA.jl/pull/135/)
- Only the relevant slices of the cached cascading matrices are loaded, instead of the whole matrices [#137](https://github.com/egavazzi/AURORA.jl/pull/137)
- Various performance improvements [#127](https://github.com/egavazzi/AURORA.jl/pull/127) [#133](https://github.com/egavazzi/AURORA.jl/pull/133) [#146](https://github.com/egavazzi/AURORA.jl/pull/146) [#151](https://github.com/egavazzi/AURORA.jl/pull/151)

## v0.7.0 - 2026-03-25
- **Breaking** Rename `animate_IeztE_3Dzoft` to `animate_Ie_in_time` [#89](https://github.com/egavazzi/AURORA.jl/pull/89)
  - Comes with a few nice improvements to the function, see PR description
- Throw an error when invalid pitch-angle limits are used as input [#90](https://github.com/egavazzi/AURORA.jl/pull/90)
- **Breaking** Add mechanism for automatic time slicing of simulations [#91](https://github.com/egavazzi/AURORA.jl/pull/91)
  - **Breaking** `calculate_e_transport()` now takes `t_total` and `dt` (in seconds) instead of `t_sampling` and `n_loop`
  - `n_loop` is now automatically calculated to keep memory usage below a configurable limit (default: 8 GB), but can still be overridden by passing it as a keyword argument
- **Breaking** Rework the input flux functions [#68](https://github.com/egavazzi/AURORA.jl/pull/68) [#109](https://github.com/egavazzi/AURORA.jl/pull/109)
  - **Breaking** Merge `Ie_top_constant()`, `Ie_top_flickering()`, and `Ie_top_Gaussian()` into a single unified `Ie_top_modulated()` function, with keyword arguments to control the energy spectrum (`:flat` or `:gaussian`) and temporal modulation (`:none`, `:sinus`, or `:square`)
  - **Breaking** `Ie_top_from_file()` has a new, simplified interface: the `n_loop` argument is removed, and the function now supports arbitrary time grids in the file (different `dt`, different length) via interpolation (`:constant`, `:linear` or `:pchip`)
  - **Breaking** `Ie_with_LET()` now takes `IeE_tot` in W/m² (instead of `Q` in eV/m²/s) as its first argument
  - **Breaking** `make_altitude_grid()` now ensures the last grid point is strictly below the requested top altitude (the grid can be one step smaller than before)
- Add possibility to save the input flux to the output directory [#103](https://github.com/egavazzi/AURORA.jl/pull/103)
  - `calculate_e_transport()` and `calculate_e_transport_steady_state()` now accept a `save_input_flux` keyword argument (default: `true`) that saves the top-boundary flux to `Ie_incoming.mat` in the output directory
- Performance improvement of animating the flux [#107](https://github.com/egavazzi/AURORA.jl/pull/107)
- Add new analysis functions to calculate phase-space density and field-aligned distribution [#108](https://github.com/egavazzi/AURORA.jl/pull/108)

## v0.6.0 - 2025-11-04
- Fix Python package installation issue with Conda [#77](https://github.com/egavazzi/AURORA.jl/pull/77)
- Make it possible to produce column excitations from steady-state results [#76](https://github.com/egavazzi/AURORA.jl/pull/76)
- Add analysis function to calculate heating rates [#73](https://github.com/egavazzi/AURORA.jl/pull/73) 
- Fix and improve the `Ie_with_LET()` function [#71](https://github.com/egavazzi/AURORA.jl/pull/68)
- Fix negative densities at very low altitudes (< 85km) [#69](https://github.com/egavazzi/AURORA.jl/pull/69)
- Increase possible maximum energy to 1 MeV (but please don't do that) [#67](https://github.com/egavazzi/AURORA.jl/pull/67)
- Improve performances many places, making simulations 5x to 15x faster to run [#64](https://github.com/egavazzi/AURORA.jl/pull/64) [#81](https://github.com/egavazzi/AURORA.jl/pull/81)
- Refactor and speed-up the cascading calculations [#72](https://github.com/egavazzi/AURORA.jl/pull/72)
- Refactor and speed-up the scattering calculations, *can change results very slightly* [#66](https://github.com/egavazzi/AURORA.jl/pull/66)


## v0.5.0 - 2025-05-01
- Faster phase functions calculations [#62](https://github.com/egavazzi/AURORA.jl/pull/62)
- Make it possible to choose a bottom altitude for the ionosphere [#58](https://github.com/egavazzi/AURORA.jl/pull/58)
- Clean the dependencies and use extensions, which reduces the precompilation times [#56](https://github.com/egavazzi/AURORA.jl/pull/56)
- Allow for saving simulation results anywhere on the system [#57](https://github.com/egavazzi/AURORA.jl/pull/57)
- Precompile some functions, leading to 10x faster simulation startup in new Julia sessions [#51](https://github.com/egavazzi/AURORA.jl/pull/51), [#52](https://github.com/egavazzi/AURORA.jl/pull/52)
- Add julia script and functions to make an animation of simulation results [#50](https://github.com/egavazzi/AURORA.jl/pull/50)
- Rewrite the analysis functions into Julia [#42](https://github.com/egavazzi/AURORA.jl/pull/42)
  - Speedup of the analysis of simulation results
  - Now load and analyze the results "slice by slice", which makes it possible to handle longer simulations
  - Emission cross-section functions translated to Julia
- Add analysis functions to the control script template 
- Remove the last Matlab dependencies [#49](https://github.com/egavazzi/AURORA.jl/pull/49)
- Performance improvement [#44](https://github.com/egavazzi/AURORA.jl/pull/44)

## v0.4.3 - 2025-01-02
- fix bug where secondary e- are not properly redistributed isotropically [#43](https://github.com/egavazzi/AURORA.jl/pull/43)
- add new docs and docstrings

## v0.4.2 - 2024-09-05
- fix bug where ionization rates have too low values due to missing secondary e- [#40](https://github.com/egavazzi/AURORA.jl/pull/40)

## v0.4.1 - 2024-07-12
- fix bug where ionization rates have too high values [#37](https://github.com/egavazzi/AURORA.jl/pull/37)

## v0.4.0 - 2024-05-22
- register the repository on [zenodo.org](https://zenodo.org/)
- add a .JuliaFormatter.toml file for the inbuilt vscode Julia extension formatter
- rewrite of the cross-section functions in Julia, which means the whole setup is now in Julia [#34](https://github.com/egavazzi/AURORA.jl/pull/34)
- add info to help debug segfault when calling Matlab [#33](https://github.com/egavazzi/AURORA.jl/pull/33)
- use 'pymsis' and 'iri2016' python packages to get msis and iri data [#30](https://github.com/egavazzi/AURORA.jl/pull/30)
- use half steps in height for A and B matrices in the CN [#28](https://github.com/egavazzi/AURORA.jl/pull/28)
- make saving simulation data safer [#27](https://github.com/egavazzi/AURORA.jl/pull/27)
- add scripts to plot I and Q in Julia [#26](https://github.com/egavazzi/AURORA.jl/pull/26)
- use a finer grid in altitude [#25](https://github.com/egavazzi/AURORA.jl/pull/25)
- add a steady state version of the transport code [#24](https://github.com/egavazzi/AURORA.jl/pull/24)

## v0.3.1 - 2023-12-18
- better calculations of dt and of the CFL factor, which yields performance improvements (see the [commit](https://github.com/egavazzi/AURORA.jl/commit/31274452819201eb28d64be530baf85cb521e291))
- fix bug https://github.com/egavazzi/AURORA.jl/issues/22

## v0.3.0 - 2023-12-14
- big performance improvements, on the order of 5x faster [#21](https://github.com/egavazzi/AURORA.jl/pull/21)
- iri data are automatically downloaded/loaded [#20](https://github.com/egavazzi/AURORA.jl/pull/20)
- nrlmsis data are automatically downloaded/loaded [#17](https://github.com/egavazzi/AURORA.jl/pull/17)
- electron densities can be calculated from Ie and can be plotted. Ionization rates can be plotted too [#15](https://github.com/egavazzi/AURORA.jl/pull/15)
- electron flux results from simulation can be downsampled in time [#16](https://github.com/egavazzi/AURORA.jl/pull/16)
- update docs about how to get started with simulations

## v0.2.0 - 2023-10-26
- the input from file function now handles non-matching time arrays [#6](https://github.com/egavazzi/AURORA.jl/pull/6)
- performance improvements of the energy degradation part [#12](https://github.com/egavazzi/AURORA.jl/pull/12)
- partial rewrite of the setup in Julia [#8](https://github.com/egavazzi/AURORA.jl/pull/8)
- the MATLAB scripts that AURORA.jl still depends on are now directly packaged with AURORA.jl [#10](https://github.com/egavazzi/AURORA.jl/pull/10)
- a bug with the calculations of beam weights and Pmu2mup matrices is fixed [#7](https://github.com/egavazzi/AURORA.jl/issues/7)
- add a proper citation file
- code is renamed to AURORA.jl
