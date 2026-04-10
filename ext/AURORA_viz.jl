module AURORA_viz

using AURORA
include("src/utilities.jl")
include("src/animations/plot.jl")
include("src/animations/animate.jl")
include("src/plot_input.jl")
include("src/plot_excitation.jl")
include("src/plot_emission.jl")
# export animate_Ie_in_time

end
