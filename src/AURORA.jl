module AURORA

# AURORA functions
include("../internal_data/data_electron/e_N2_cross_sections.jl")
export e_N2elastic, e_N2rot0_2, e_N2rot0_4, e_N2rot0_6, e_N2rot0_8, e_N2vib0_1, e_N2vib0_2,
    e_N2vib0_3, e_N2vib0_4, e_N2vib0_5, e_N2vib0_6, e_N2vib0_7, e_N2a3sup, e_N2b3pg,
    e_N2w3du, e_N2bp3sum, e_N2ap1sum, e_N2w1du, e_N2e3sgp, e_N2ab1sgp, e_N2a1pg, e_N2c3pu,
    e_N2bp1sup, e_N2cp1sup, e_N2cp3pu, e_N2d3sup, e_N2f3pu, e_N2g3pu, e_N2M1M2, e_N2o1pu,
    e_N2dissociation, e_N2ionx2sgp, e_N2iona2pu, e_N2ionb2sup, e_N2dion, e_N2ddion
include("../internal_data/data_electron/e_O2_cross_sections.jl")
export e_O2elastic, e_O2_OO3S, e_O2_9p97, e_O2_8p4, e_O2_6, e_O2_4p5, e_O2b1Sgp, e_O2a1Dg,
    e_O2vib, e_O2ionx2pg, e_O2iona4pu, e_O2ion16p9, e_O2ionb4sgm, e_O2dion, e_O2ddion
include("../internal_data/data_electron/e_O_cross_sections.jl")
export e_Oelastic, e_O1D, e_O1S, e_O3s5S0, e_O3s3S0, e_O3p5P, e_O3sp3D0, e_O3p3P, e_Oion4S0,
    e_Oion2D0, e_Oion2P0, e_Oionion

include("setup.jl")
include("nrlmsis.jl")
include("iri.jl")
include("input_flux.jl")
include("phase_functions.jl")
include("utilitaries.jl")
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
       f_smooth_transition, rename_if_exists
export d2M, Crank_Nicolson, Crank_Nicolson_Optimized
export cascading_N2, cascading_O2, cascading_O
export update_Q!, update_Q_turbo!
export rotating, load_scattering_matrices

include("main.jl")
export calculate_e_transport

include("steady_state.jl")
export calculate_e_transport_steady_state

include("analysis.jl")
export make_density_file, downsampling_fluxes



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
