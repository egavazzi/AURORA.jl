module Aurora

using MATLAB

include("setup.jl")
include("matrix_building.jl")
include("phase_functions.jl")
include("utilitaries.jl")

export setup
export phase_fcn_N2, phase_fcn_O2, phase_fcn_O, convert_phase_fcn_to_3D
export beams2beams, make_A, make_B
export v_of_E

end
