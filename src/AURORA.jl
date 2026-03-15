module AURORA

# AURORA functions


include("grids/abstract_grid.jl")
include("grids/altitude_grid.jl")
include("grids/energy_grid.jl")
include("grids/pitch_angle_grid.jl")
export AbstractGrid, AltitudeGrid, EnergyGrid, PitchAngleGrid

include("setup.jl")
export setup, make_altitude_grid, make_energy_grid

include("ionosphere/ionosphere.jl")
export Ionosphere, n_neutrals
include("ionosphere/iri/iri.jl")
include("ionosphere/msis/msis.jl")
export find_msis_file, find_nrlmsis_file
export find_iri_file

include("input/input_flux.jl")
export Ie_top_from_file, Ie_top_modulated, Ie_with_LET

include("utilities.jl")
export v_of_E, CFL_criteria, mu_avg, beam_weight,
        make_savedir,
        rename_if_exists,
        find_Ietop_file

include("solvers/transport_matrices.jl")
include("solvers/matrix_building.jl")
include("solvers/crank_nicolson.jl")
include("solvers/crank_nicolson_optimized.jl")

include("physics/cross_sections/e_N2_cross_sections.jl")
include("physics/cross_sections/e_O2_cross_sections.jl")
include("physics/cross_sections/e_O_cross_sections.jl")
include("physics/cross_sections/emission_cross_sections.jl")
include("physics/cross_sections/cross_sections.jl")
export CrossSectionData
include("physics/cascading.jl")
include("physics/energy_degradation.jl")
include("physics/scattering.jl")
export ScatteringData
include("physics/phase_functions.jl")
export phase_fcn_N2, phase_fcn_O2, phase_fcn_O, convert_phase_fcn_to_3D

include("main.jl")
export calculate_e_transport

include("solvers/steady_state.jl")
include("solvers/steady_state_optimized.jl")
export calculate_e_transport_steady_state

include("analysis.jl")
export make_density_file, downsampling_fluxes, make_volume_excitation_file,
    make_column_excitation_file, make_Ie_top_file, make_current_file, make_heating_rate_file

# Define and export functions to be extented by the AURORA_viz module
"""
    animate_Ie_in_time(directory_to_process; angles_to_plot=nothing, colorrange=nothing, ...)

Plot a heatmap of Ie over height and energy, and animate it in time. It will load the
result files one by one. The animation will be saved as a .mp4 file under the
`directory_to_process`.

# Example
```julia-repl
julia> directory_to_process = "Visions2/Alfven_475s";

# Using defaults for angles and colorrange:
julia> animate_Ie_in_time(directory_to_process)

# Or with custom angles and colorrange:
julia> angles_to_plot = [(180, 170)  (170, 150)  (150, 120)  (120, 100)  (100, 90);   # DOWN
                         (0, 10)     (10, 30)    (30, 60)    (60, 80)    (80, 90)];   # UP
julia> animate_Ie_in_time(directory_to_process; angles_to_plot, colorrange=(1e5, 1e9), plot_Ietop=true)

# Using nothing for empty panels:
julia> angles_to_plot = [(180, 90)  nothing;
                         (0, 45)    (45, 90)];
julia> animate_Ie_in_time(directory_to_process; angles_to_plot)
```

The `angles_to_plot` is a matrix of tuples, where each tuple defines a pitch-angle range
from 0° to 180° (where 180° is field-aligned down and 0° is field-aligned up). A panel
will be created for each matrix element at the corresponding row/column position.
Angles > 90° are labeled as "DOWN", angles < 90° as "UP". Use `nothing` for empty panels.

The limits of `angles_to_plot` need to match existing limits of the beams
used in the simulation. E.g. if `θ_lims = 180:-10:0` was used in the simulation, `(150, 120)`
will be fine as 150° and 120° exist as limits, but `(155, 120)` will not as 155° does not
exist as a limit.

# Arguments
- `directory_to_process`: directory containing the simulation results (absolute or relative path).

# Keyword Arguments
- `angles_to_plot = nothing`: limits of the angles to plot as a matrix of tuples with angles
                              in range 0-180°. Use `nothing` for empty panels. If the whole
                              argument is `nothing`, uses the θ_lims grid from the simulation
                              with down-flux on the first row and up-flux on the second row.
- `colorrange = nothing`: limits for the colormap/colorbar as a tuple (min, max). If `nothing`,
                          automatically computed as (max_value / 1e4, max_value) spanning
                          4 orders of magnitude.
- `save_to_file = true`: if true, saves the animation to a .mp4 file in the data directory.
- `plot_Ietop = false`: if true, also plots the precipitating Ie at the top of the
                        ionosphere by loading it from the file `Ie_top.mat`.
- `Ietop_angle_cone = [170, 180]`: angle cone (in degrees) for the precipitating Ie
                        to plot.
"""
function animate_Ie_in_time end
export animate_Ie_in_time

# Precompile selected functions
include("precompiles.jl")

end
