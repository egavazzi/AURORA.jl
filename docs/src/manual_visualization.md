# Visualization

AURORA.jl provides functions to visualize and animate simulation results. 

The main animation
function is `animate_Ie_in_time`.

## Requirements

The plotting and animation functions require a [Makie](https://github.com/MakieOrg/Makie.jl) backend.

We recommend:
- **GLMakie** - if you have a GPU, is fast and interactive
- **CairoMakie** - if you work on a server without a GPU or want high quality figures (e.g. for publication). It is however non-interactive and slower


!!! note "Installing Makie"
    It can be a good idea to install GLMakie or CairoMakie in your global Julia environment (`@v1.XX`)

To use the visualization functions, call a backend before or after calling AURORA
```julia
# For example
using AURORA
using GLMakie

# Do some visualizaton
```


## Functions

```@docs
animate_Ie_in_time
```