#=
Functions to animate Ie(z, t, E).

So far we have
- animate_IeztE_3Dzoft: Ie as a heatmap over height and energy, animation in time.


=#
using GLMakie
using CairoMakie
using MAT
GLMakie.activate!()

# First a function to load Ie
filename = "/run/user/1000/gvfs/sftp:host=revontuli.uit.no/mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_536s_longer_fixed-secondaries/IeFlickering-01.mat"
filename = "data/Visions2/IeFlickering-01.mat"
@time data = matread(filename)
@time Ie_raw = data["Ie_ztE"] # size of [n_mu x nz, nt, nE]

μ_scatterings = data["mu_scatterings"]
μ_lims = data["mu_lims"]
t_run = data["t_run"]
E = data["E"]
h_atm = data["h_atm"]

## Then some function to restructure the matrix from [n_mu x nz, nt, nE] to [n_mu, nz, nt, nE]
n_z = length(h_atm)
n_μ = length(μ_lims) - 1
n_t = length(t_run)
n_E = length(E)
i_μ = 1

function restructure_Ie(Ie_raw, μ_lims, h_atm, t_run, E)
    n_μ = length(μ_lims) - 1
    n_z = length(h_atm)
    n_t = length(t_run)
    n_E = length(E)
    Ie_restructured = zeros(n_μ, n_z, n_t, n_E);
    for i_E in 1:n_E
        for i_t in 1:n_t
            for i_z in 1:n_z
                for i_μ in 1:n_μ
                    Ie_restructured[i_μ, i_z, i_t, i_E] = Ie_raw[i_z + (i_μ - 1) * n_z, i_t, i_E]
                end
            end
        end
    end
    return Ie_restructured # size [n_mu, nz, nt, nE]
end

@time Ie = restructure_Ie(Ie_raw, μ_lims, h_atm, t_run, E); # size [n_mu, nz, nt, nE]

Ie_timeslice = Observable(Ie[:, :, 1, :])
time = Observable("$(t_run[1]) s")

##
#=
This is the function that creates the plot from the data observable. Then when we update the
data, the plot is automatically updated.

It should take Ie as a function of pitch-angle, height and energy as input. But NOT as a
function of time. This because we want to update the time OUTSIDE of the function.

TODO: Add a function before that one that remodel Ie with the pitch-angle streams that we
want to plot. Then we need to change the function a bit to have one axis per pitch-angle
stream.
=#
function make_IeztE_3Dzoft_plot(Ie_timeslice::Observable{Array{Float64, 3}}, time::Observable{String}, limits)
    # Slice the input Ie into its different pitch-angle components
    Ie_1 = @lift($Ie_timeslice[1, :, :]')
    Ie_2 = @lift($Ie_timeslice[2, :, :]')
    Ie_3 = @lift($Ie_timeslice[3, :, :]')
    Ie_4 = @lift($Ie_timeslice[4, :, :]')

    # Plot (redo with loop and proper subplots down and up)
    f = Figure(size = (600, 450))
    # ax1 = Axis(f[1, 1]; xscale = log10, title = time)
    ax1 = Axis(f[1, 1]; xscale = log10)
    hm1 = heatmap!(E, h_atm / 1e3, Ie_1, colorrange = limits, colorscale = log10)

    ax2 = Axis(f[2, 1]; xscale = log10)
    hm2 = heatmap!(E, h_atm / 1e3, Ie_2, colorrange = limits, colorscale = log10)

    ax3 = Axis(f[3, 1]; xscale = log10)
    hm3 = heatmap!(E, h_atm / 1e3, Ie_3, colorrange = limits, colorscale = log10)

    ax4 = Axis(f[4, 1]; xscale = log10)
    hm4 = heatmap!(E, h_atm / 1e3, Ie_4, colorrange = limits, colorscale = log10)

    cb = Colorbar(f[1, 2], hm1)
    Label(f[0, 1], time; tellwidth = false, fontsize=18)

    return f
end

## Call the function

global_min, global_max = extrema(Ie)
global_min = 1e0

f = make_IeztE_3Dzoft_plot(Ie_timeslice, time, (global_min, global_max))
display(f)

## Test the animation
# @time for i_t in 1:10
#     Ie_timeslice[] .= Ie[:, :, i_t, :]
#     notify(Ie_timeslice)
#     sleep(0.01)
# end

## Record the animation to a mp4
# with glmakie(~12s for n_t=51)
@time record(f, "time_animation1_v2.mp4", 1:n_t; framerate=10, backend=GLMakie) do i_t
    Ie_timeslice[] .= Ie[:, :, i_t, :]
    time[] = "$(t_run[i_t]) s"
    notify(Ie_timeslice)
end
# with cairomakie (~110s for n_t=51)
@time record(f, "time_animation2_v2.mp4", 1:n_t; framerate=10, backend=CairoMakie) do i_t
    Ie_timeslice[] .= Ie[:, :, i_t, :]
    time[] = "$(t_run[i_t]) s"
    notify(Ie_timeslice)
end


##
# Then a function to animate Ie in time
function animate_IeztE_3Dzoft(t, h_atm, E, dE, Ie)
    # Ie (in #e⁻/m²/s/eV/ster)
    f = Figure()
    ax = Axis
end
