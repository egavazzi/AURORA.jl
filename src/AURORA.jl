module AURORA

# AURORA functions
include("../internal_data/data_electron/e_N2_cross_sections.jl")
include("../internal_data/data_electron/e_O2_cross_sections.jl")
include("../internal_data/data_electron/e_O_cross_sections.jl")
include("../internal_data/data_electron/emission_cross_sections.jl")

include("setup.jl")
include("nrlmsis.jl")
include("iri.jl")
include("input_flux.jl")
include("phase_functions.jl")
include("utilities.jl")
include("matrix_building.jl")
include("crank_nicolson.jl")
include("cascading.jl")
include("energy_degradation.jl")
include("scattering.jl")
export setup, make_altitude_grid, make_energy_grid, make_scattering_matrices
export find_nrlmsis_file
export find_iri_file
export Ie_top_from_file, Ie_top_flickering, Ie_top_constant
export phase_fcn_N2, phase_fcn_O2, phase_fcn_O, convert_phase_fcn_to_3D
export loss_to_thermal_electrons, beams2beams, make_A, make_B, make_D
export v_of_E, CFL_criteria, mu_avg, beam_weight, save_parameters, save_results,
       f_smooth_transition, rename_if_exists, find_Ietop_file, make_savedir
export d2M, Crank_Nicolson
export cascading_N2, cascading_O2, cascading_O
export update_Q!, update_Q_turbo!
export rotating, load_scattering_matrices

include("main.jl")
export calculate_e_transport

include("steady_state.jl")
export calculate_e_transport_steady_state

include("analysis.jl")
export make_density_file, downsampling_fluxes, make_volume_excitation_file,
    make_column_excitation_file, make_Ie_top_file, make_current_file

# Define and export functions to be extented by the AURORA_viz module
function animate_IeztE_3Dzoft end
export animate_IeztE_3Dzoft

# Precompile selected functions
include("precompiles.jl")

end
