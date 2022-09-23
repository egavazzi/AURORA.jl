module Aurora

include("setup.jl")
include("input_flux.jl")
include("phase_functions.jl")
include("utilitaries.jl")
include("matrix_building.jl")
include("crank_nicolson.jl")
include("cascading.jl")
include("energy_degradation.jl")

export setup
export Ie_top_from_file
export phase_fcn_N2, phase_fcn_O2, phase_fcn_O, convert_phase_fcn_to_3D
export loss_to_thermal_electrons, beams2beams, make_A, make_B, make_D
export v_of_E
export Crank_Nicolson, d2M, Crank_Nicolson_Optimized
export cascading_N2, cascading_O2, cascading_O
export update_Q!

end