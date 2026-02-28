# Visualization

AURORA.jl provides functions to visualize and animate simulation results. 

Currently, the only animation function available is `animate_Ie_in_time`.

## Requirements

The plotting and animation functions require a [Makie](https://github.com/MakieOrg/Makie.jl) backend.

We recommend:
- **GLMakie** - if you have a GPU. It is fast and interactive
- **CairoMakie** - if you work on a server without a GPU or want high quality figures (e.g. for publication). It is however non-interactive and slower than GLMakie


!!! note "Installing Makie"
    You need to install GLMakie or CairoMakie yourself, for example 
    in your global Julia environment (`@v1.XX`)

To use the visualization functions, import a backend before or after importing AURORA
```julia
# For example
using AURORA
using GLMakie

# Do some visualization
```


## Functions

```@docs; canonical=false
animate_Ie_in_time
```