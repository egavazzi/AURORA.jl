#=
Functions to animate Ie(z, t, E).

So far we have
- animate_Ie_in_time: Ie as a heatmap over height and energy, animation in time.


=#
using AURORA

#=
IF YOU HAVE A GPU, you can use GLMakie which will be much faster for plotting the animation.
Note that GLMakie needs to be installed in your global julia environment (@v1.XX).

IF YOU DON'T HAVE A GPU (for example you are running the code on a server), you need to use
CairoMakie. Comment out the lines with GLMakie and uncomment the lines with CairoMakie.
Note that this will be quite slow.
=#
using GLMakie
GLMakie.activate!()
# using CairoMakie
# CairoMakie.activate!()


## Calling the animate function
directory_to_process = "pitch-angle-tests/180:-45:0"

# animate_Ie_in_time(directory_to_process; angles_to_plot, colorrange = (1e5, 1e9), plot_Ietop = true)
animate_Ie_in_time(directory_to_process; save_to_file = true)
