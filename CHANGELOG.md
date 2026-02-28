# Changelog

- **Breaking** Rename `animate_IeztE_3Dzoft` to `animate_Ie_in_time` [#89](https://github.com/egavazzi/AURORA.jl/pull/89)
  - Comes with a few nice improvements to the function, see PR description
- Throw an error when invalid pitch-angle limits are used as input [#90](https://github.com/egavazzi/AURORA.jl/pull/90)
- **Breaking** Add mechanism for automatic time slicing of simulations [#91](https://github.com/egavazzi/AURORA.jl/pull/91)
  - **Breaking** `calculate_e_transport()` now takes `t_total` and `dt` (in seconds) instead of `t_sampling` and `n_loop`
  - `n_loop` is now automatically calculated to keep memory usage below a configurable limit (default: 8 GB), but can still be overridden by passing it as a keyword argument
- **Breaking** Rework the input flux functions [#68](https://github.com/egavazzi/AURORA.jl/pull/68)
  - **Breaking** Merge `Ie_top_constant()`, `Ie_top_flickering()`, and `Ie_top_Gaussian()` into a single unified `Ie_top_modulated()` function, with keyword arguments to control the energy spectrum (`:flat` or `:gaussian`) and temporal modulation (`:none`, `:sinus`, or `:square`)
  - **Breaking** `Ie_top_from_file()` has a new, simplified interface: the `n_loop` argument is removed, and the function now supports arbitrary time grids in the file (different `dt`, different length) via interpolation (`:constant` or `:linear`)
  - **Breaking** `Ie_with_LET()` now takes `IeE_tot` in W/m² (instead of `Q` in eV/m²/s) as its first argument
  - **Breaking** `make_altitude_grid()` now ensures the last grid point is strictly below the requested top altitude (the grid can be one step smaller than before)

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
