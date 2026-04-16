module AURORA

include("grids/abstract_grid.jl")
include("grids/altitude_grid.jl")
include("grids/energy_grid.jl")
include("grids/pitch_angle_grid.jl")
export AbstractGrid, AltitudeGrid, EnergyGrid, PitchAngleGrid

include("ionosphere/ionosphere.jl")
export Ionosphere, n_neutrals
include("ionosphere/iri/iri.jl")
include("ionosphere/msis/msis.jl")
export find_msis_file, find_nrlmsis_file
export find_iri_file

include("physics/cross_sections/e_N2_cross_sections.jl")
include("physics/cross_sections/e_O2_cross_sections.jl")
include("physics/cross_sections/e_O_cross_sections.jl")
include("physics/cross_sections/emission_cross_sections.jl")
include("physics/cross_sections/cross_sections.jl")
export CrossSectionData

include("physics/scattering.jl")
export ScatteringData
include("physics/phase_functions.jl")
export phase_fcn_N2, phase_fcn_O2, phase_fcn_O, convert_phase_fcn_to_3D

include("model.jl")
export AuroraModel, make_altitude_grid, make_energy_grid

include("input/spectra.jl")
include("input/modulations.jl")
include("input/input_flux.jl")
export AbstractSpectrum, FlatSpectrum, GaussianSpectrum, MaxwellianSpectrum, FileSpectrum
export AbstractModulation, ConstantModulation, SinusoidalFlickering, SquareFlickering, SmoothOnset
export InputFlux, evaluate_spectrum, apply_modulation, compute_flux
export Ie_top_from_file


include("solvers/transport_matrices.jl")
include("solvers/matrix_building.jl")

include("physics/cascading.jl")
include("physics/energy_degradation.jl")

include("solvers/crank_nicolson.jl")
include("solvers/crank_nicolson_optimized.jl")
include("solvers/steady_state.jl")
include("solvers/steady_state_optimized.jl")

include("simulation/cache.jl")
include("simulation/types.jl")
export AuroraSimulation, ResolvedTimeGrid
include("simulation/initialize.jl")
export initialize!
include("simulation/run.jl")
export run!

include("utilities.jl")
export v_of_E, mu_avg, beam_weight, make_savedir

include("analysis/analysis_types.jl")
export VolumeExcitationResult, ColumnExcitationResult, IeTopResult,
       load_volume_excitation, load_column_excitation, load_input
include("analysis/psd.jl")
include("analysis/emissions.jl")
include("analysis/fluxes.jl")
include("analysis/heating.jl")
export make_volume_excitation_file, make_column_excitation_file,
       downsampling_fluxes, make_Ie_top_file, make_current_file,
       make_heating_rate_file, make_psd_file

# Define and export functions to be extended by the AURORA_viz module
"""
    plot_input(sim::AuroraSimulation)

Plot the input electron flux for a simulation.

For **time-dependent** simulations, produces a heatmap of flux vs energy and time for each
active beam.  For **steady-state** simulations, produces a line plot of flux vs energy.

Requires a Makie backend (e.g. `using CairoMakie` or `using GLMakie`).

# Examples
```julia
using CairoMakie
sim = AuroraSimulation(model, flux, 0.5, 0.001, savedir)
fig = plot_input(sim)
```
"""
function plot_input end
"""
    plot_input!(ax, data::IeTopResult;
                beams=1, colorrange=nothing, colormap=:inferno, kwargs...)

Plot the input electron flux stored in an [`IeTopResult`](@ref) onto an existing `Axis`.

Requires a Makie backend (e.g. `using CairoMakie` or `using GLMakie`).

# Keyword Arguments
- `beams = 1`: pitch-angle beam index, or vector of indices, to include. Selected beams are
  summed and normalised by their combined solid angle, then converted to differential energy
  flux (#e⁻/m²/s/eV/ster).
- `colorrange = nothing`: `(min, max)` limits for the colorscale. Defaults to
  `(max_value / 1e4, max_value)` spanning 4 orders of magnitude.
- `colormap = :inferno`: Makie colormap for the heatmap.
- `kwargs...`: additional keyword arguments are forwarded to the underlying `heatmap!` call.

"""
function plot_input! end
export plot_input, plot_input!

"""
    animate_Ie_in_time(directory_to_process;
                       angles_to_plot=nothing, colorrange=nothing, save_to_file=true,
                       plot_input=false, input_angle_cone=[170, 180], dt_steps=1,
                       framerate=30)

Plot a heatmap of Ie over height and energy, and animate it in time. The animation is saved
as a .mp4 file under the `directory_to_process` if `save_to_file = true`.

Requires a Makie backend (e.g. `using CairoMakie` or `using GLMakie`).

# Arguments
- `directory_to_process`: directory containing the simulation results (absolute or relative path).

# Keyword Arguments
- `angles_to_plot = nothing`: limits of the angles to plot as a matrix of tuples with angles
                              in range 0-180°. Use `nothing` for empty panels. If the whole
                              argument is `nothing`, uses the simulation grid
                              with down-flux on the first row and up-flux on the second row.
- `colorrange = nothing`: limits for the colormap/colorbar as a tuple (min, max). If `nothing`,
                          automatically computed as (max_value / 1e4, max_value) spanning
                          4 orders of magnitude.
- `save_to_file = true`: if `true`, saves the animation to an `animation.mp4` file (or
                         `animation_with_precipitation.mp4` when `plot_input = true`) inside
                         `directory_to_process`.
- `plot_input = false`: if `true`, also plots the precipitating electron flux at the top of
                        the ionosphere by loading it from the `Ie_incoming_*.mat` file.
- `input_angle_cone = [170, 180]`: pitch-angle cone (degrees) used to select and sum beams
                                   for the precipitation overlay. Only used when `plot_input = true`.
- `dt_steps = 1`: plot one frame every `dt_steps` timesteps. Increase to speed up rendering.
- `framerate = 30`: framerate of the animation in frames per second.

# Example
```julia-repl
# Use default angles and colorrange:
julia> animate_Ie_in_time(directory_to_process)

# Use custom angles and colorrange:
julia> angles_to_plot = [(180, 170)  (170, 150)  (150, 120)  (120, 100)  (100, 90);   # DOWN
                         (0, 10)     (10, 30)    (30, 60)    (60, 80)    (80, 90)];   # UP
julia> animate_Ie_in_time(directory_to_process; angles_to_plot, colorrange=(1e5, 1e9), plot_input=true)

# Use `nothing` for an empty panels:
julia> angles_to_plot = [(180, 90)  nothing;
                         (0, 45)    (45, 90)];
julia> animate_Ie_in_time(directory_to_process; angles_to_plot)
```

# Notes
The `angles_to_plot` is a matrix of tuples, where each tuple defines a pitch-angle range
from 0° to 180° (where 180° is field-aligned down and 0° is field-aligned up). A panel
will be created for each matrix element at the corresponding row/column position.
Angles > 90° are labeled as "DOWN", angles < 90° as "UP". Use `nothing` for empty panels.

The limits of `angles_to_plot` need to match existing limits of the beams
used in the simulation. E.g. if `θ_lims = 180:-10:0` was used in the simulation, `(150, 120)`
will be fine as 150° and 120° exist as limits, but `(155, 120)` will not as 155° does not
exist as a limit.
"""
function animate_Ie_in_time end
export animate_Ie_in_time

"""
    plot_excitation!(ax, data::VolumeExcitationResult;
                     field=:total, time_index=nothing, colorrange=nothing,
                     colormap=:viridis, show_contours=false, show_height_of_max=true,
                     kwargs...)

Plot a volume excitation or ionization rate onto an existing `Axis`.

Requires a Makie backend (e.g. `using CairoMakie` or `using GLMakie`).


- If `data` contains one time step, or if `time_index` is given, plots an altitude
  profile (returns a `Lines` plot).
- Otherwise, plots a time-altitude heatmap (returns a `Heatmap`).

# Keyword Arguments
- `field = :total`: which rate to plot. Can be any field of [`VolumeExcitationResult`](@ref):
  `:Q4278`, `:Q6730`, `:Q7774`, `:Q8446`, `:QO1D`, `:QO1S`, `:QOi`, `:QO2i`, `:QN2i`,
  or `:total` which sums the three ionization rates (`:QN2i + :QO2i + :QOi`).
- `time_index = nothing`: if given, plot a single altitude profile at that time index.
  If `nothing` and the data has more than one time step, a heatmap is plotted.
- `colorrange = nothing`: `(min, max)` limits for the colorscale. Defaults to
  `(max_value / 1e4, max_value)` spanning 4 orders of magnitude.
- `colormap = :viridis`: Makie colormap for the heatmap.
- `show_contours = false`: if `true`, overlays black contour lines on the heatmap (heatmap mode only).
- `show_height_of_max = true`: if `true`, overlays a dashed red line at the altitude of peak
  rate at each time step (heatmap mode only).
- `kwargs`: additional keyword arguments are forwarded to the underlying `heatmap!` or `lines!` call.
"""
function plot_excitation! end
export plot_excitation!

"""
    plot_column_excitation!(ax, data::ColumnExcitationResult;
                            wavelengths=[:I_4278, :I_6730, :I_7774, :I_8446, :I_O1D, :I_O1S],
                            kwargs...)

Plot column-integrated emission intensities onto an existing `Axis`.

Requires a Makie backend (e.g. `using CairoMakie` or `using GLMakie`).

# Keyword Arguments
- `wavelengths`: vector of symbols selecting which emission lines to plot. Available values
  are `:I_4278`, `:I_6730`, `:I_7774`, `:I_8446`, `:I_O1D`, `:I_O1S`. Defaults to all six.
  O(¹D) and O(¹S) lines are plotted with a dashed linestyle.
- `kwargs`: additional keyword arguments are forwarded to each underlying `lines!` call.
"""
function plot_column_excitation! end
export plot_column_excitation!

"""
    plot_model(model::AuroraModel; panels=[:all])

Plot the model setup: atmosphere, energy grid, cross-sections, phase functions, and
scattering data.

Returns a `Dict{Symbol, Figure}` with one figure per panel. Use the `panels` keyword to
select which panels to plot.

Requires a Makie backend (e.g. `using CairoMakie` or `using GLMakie`).

# Available panels
`:atmosphere`, `:energy_levels`, `:energy_grid`, `:cross_sections`, `:phase_functions`,
`:scattering`, `:beam_weights`, or `:all` (default) to plot everything.

# Examples
```julia
using CairoMakie
figs = plot_model(model)
figs = plot_model(model; panels=[:atmosphere, :cross_sections])
```
"""
function plot_model end
export plot_model

# Precompile selected functions
include("precompiles.jl")

end
