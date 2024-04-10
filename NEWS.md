# News

- use 'pymsis' and 'iri2016' python packages to get msis and iri data [#30](https://github.com/egavazzi/AURORA.jl/pull/30)
- add a .JuliaFormatter.toml file for the inbuilt vscode Julia extension formatter
- use half steps in height for A and B matrices in the CN [#28](https://github.com/egavazzi/AURORA.jl/pull/28)
- make saving simulation data safer [#27](https://github.com/egavazzi/AURORA.jl/pull/27)
- add scripts to plot I and Q in Julia [#26](https://github.com/egavazzi/AURORA.jl/pull/26)
- use a finer grid in altitude [#25](https://github.com/egavazzi/AURORA.jl/pull/25)
- add a steady state version of the transport code [#24](https://github.com/egavazzi/AURORA.jl/pull/24)

## v0.3.1
- better calculations of dt and of the CFL factor, which yields performance improvements (see the [commit](https://github.com/egavazzi/AURORA.jl/commit/31274452819201eb28d64be530baf85cb521e291))
- fix bug https://github.com/egavazzi/AURORA.jl/issues/22

## v0.3.0
- big performance improvements, on the order of 5x faster [#21](https://github.com/egavazzi/AURORA.jl/pull/21)
- iri data are automatically downloaded/loaded [#20](https://github.com/egavazzi/AURORA.jl/pull/20)
- nrlmsis data are automatically downloaded/loaded [#17](https://github.com/egavazzi/AURORA.jl/pull/17)
- electron densities can be calculated from Ie and can be plotted. Ionization rates can be plotted too [#15](https://github.com/egavazzi/AURORA.jl/pull/15)
- electron flux results from simulation can be downsampled in time [#16](https://github.com/egavazzi/AURORA.jl/pull/16)
- update docs about how to get started with simulations

## v0.2.0
- the input from file function now handles non-matching time arrays [#6](https://github.com/egavazzi/AURORA.jl/pull/6)
- performance improvements of the energy degradation part [#12](https://github.com/egavazzi/AURORA.jl/pull/12)
- partial rewrite of the setup in Julia [#8](https://github.com/egavazzi/AURORA.jl/pull/8)
- the MATLAB scripts that AURORA.jl still depends on are now directly packaged with AURORA.jl [#10](https://github.com/egavazzi/AURORA.jl/pull/10)
- a bug with the calculations of beam weights and Pmu2mup matrices is fixed [#7](https://github.com/egavazzi/AURORA.jl/issues/7)
- add a proper citation file
- code is renamed to AURORA.jl
