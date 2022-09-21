module Aurora

include("setup.jl")
include("matrix_building.jl")
include("phase_functions.jl")
include("utilitaries.jl")
include("crank_nicolson.jl")
include("cascading.jl")
include("energy_degradation.jl")

export setup
export phase_fcn_N2, phase_fcn_O2, phase_fcn_O, convert_phase_fcn_to_3D
export loss_to_thermal_electrons, beams2beams, make_A, make_B, make_D
export v_of_E, cascading_N2, cascading_O2, cascading_O
export Crank_Nicolson, d2M
export update_Q!, update_Q2!

end
