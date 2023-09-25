module AURORA

include("../internal_data/data_electron/e_N2_cross_sections.jl")
export e_N2elastic, e_N2rot0_2

include("setup.jl")
include("input_flux.jl")
include("phase_functions.jl")
include("utilitaries.jl")
include("matrix_building.jl")
include("crank_nicolson.jl")
include("cascading.jl")
include("energy_degradation.jl")
include("scattering.jl")
export setup, setup_new
export Ie_top_from_file, Ie_top_flickering, Ie_top_constant
export phase_fcn_N2, phase_fcn_O2, phase_fcn_O, convert_phase_fcn_to_3D
export loss_to_thermal_electrons, beams2beams, make_A, make_B, make_D
export v_of_E, CFL_criteria, mu_avg, beam_weight, save_parameters, save_results, f_smooth_transition
export d2M, Crank_Nicolson, Crank_Nicolson_Optimized
export cascading_N2, cascading_O2, cascading_O
export update_Q!
export rotating, load_scattering_matrices

include("main.jl")
export calculate_e_transport






# MI_coupling functions
include("../MI_coupling/src/utilitaries.jl")
include("../MI_coupling/src/ketchup_conversion.jl")
include("../MI_coupling/src/conversions.jl")
include("../MI_coupling/src/make_Ie_from_ketchup.jl")
include("../MI_coupling/src/make_f_from_AURORA.jl")

export load_fzvzmu_parallel, load_fzvzmu_serial, load_fzvzmuIB_serial, load_Bfield
export findnearestindex
export read_input
export convert_fzvzmu_to_Ie, convert_M_to_I, convert_Ie_to_fzvzmu, convert_I_to_M
export make_Ie_from_ketchup
export make_f_from_AURORA

end
