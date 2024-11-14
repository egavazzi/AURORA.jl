#=
Functions to animate Ie(z, t, E).

So far we have
- animate_IeztE_3Dzoft: Ie as a heatmap over height and energy, animation in time.


=#
using GLMakie
using CairoMakie
using MAT
using AURORA
GLMakie.activate!()


# Function to restructure the matrix from 3D [n_mu x nz, nt, nE] to 4D [n_mu, nz, nt, nE]
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


"""
Function that merges the streams/beams of Ie with θ_lims of the simulation to fit the
new_θ_lims of interest. It is useful for plotting for example.

For example, if we have θ_lims
    [180 160 140 120 100 90 80 60 40 20 0] # simulation
and we want to plot with new_θ_lims
    [180 160 120 100 90 80 60 40 20 0] # to plot
the function will merge the streams with index 2 (160-140°) and 3 (140-120°) together.

# Arguments
- `Ie`: electron flux from the simulation, with pitch-angle limits θ_lims. Of shape [n_μ, n_z, n_t, n_E].
- `θ_lims`: pitch-angle limits used of the simulation.
- `new_θ_lims`: new pitch-angle limits. Given as an array of tuples, for example:
                ```julia-repl
                new_θ_lims = [(0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90);  # DOWN
                              (0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90)]  # UP
                ```

# Output
- `Ie_plot`: electron flux restructured to fit the new pitch-angle limits. Of shape [n_μ_new, n_z, n_t, n_E], where
n_μ_new is the number of elements in new_θ_lims. The first dimension (μ) of Ie_plot is sorted such that the
pitch-angle stream indices go along the first row of new_θ_lims, and then the second row.
In our example with new_θ_lims, that's [1 2 3 4 5;
                                        6 7 8 9 10]
"""
function restructure_streams_of_Ie(Ie, θ_lims, new_θ_lims)

    # Initialize the new Ie that will contain the restructured streams
    n_μ_new = length(new_θ_lims)
    n_z = size(Ie, 2)
    n_t = size(Ie, 3)
    n_E = size(Ie, 4)
    Ie_plot = zeros(n_μ_new, n_z, n_t, n_E)

    # Modify the values of the down-angles so that field aligned is 180°.
    # For example new_θ_lims would go from something like
    #   2×5 Matrix{Tuple{Int64, Int64}}:
    #   (0, 10)  (10, 30)  (30, 60)  (60, 80)  (80, 90) # DOWN
    #   (0, 10)  (10, 30)  (30, 60)  (60, 80)  (80, 90) # UP
    # to
    #   2×5 Matrix{Tuple{Int64, Int64}}:
    #   (180, 170)  (170, 150)  (150, 120)  (120, 100)  (100, 90) # DOWN
    #   (0, 10)     (10, 30)    (30, 60)    (60, 80)    (80, 90)  # UP
    new_θ_lims_temp = copy(new_θ_lims)
    new_θ_lims_temp[1, :] = map(x -> 180 .- x, new_θ_lims_temp[1, :] )


    # Restructure to [n_μ_new, n_z, n_t, n_E] --- METHOD 1
    n_row = size(new_θ_lims, 1)
    n_col = size(new_θ_lims, 2)
    @views for i in axes(new_θ_lims, 1)
        for j in axes(new_θ_lims, 2)
            idx = (i - 1) * n_col  + j # goes along first row, then second row
            # for i_μ_old in axes(Ie, 1)
            #     @views if minimum(new_θ_lims_temp[i, j]) .<= acosd.(mu_avg(θ_lims))[i_μ_old] .<= maximum(new_θ_lims_temp[i, j])
            #         Ie_plot[idx, :, :, :] .+= Ie[i_μ_old, :, :, :]
            #     end
            # end
            idx_θ = minimum(new_θ_lims_temp[i, j]) .<= acosd.(mu_avg(θ_lims)) .<= maximum(new_θ_lims_temp[i, j])
            Ie_plot[idx, :, :, :] .= dropdims(sum(Ie[idx_θ, :, :, :]; dims = 1); dims = 1)
        end
    end

    # Restructure to [n_μ_new, n_z, n_t, n_E] --- METHOD 2
    # Restructure `new_θ_lims` from a 2D array to a 1D vector, by concatenating the rows.
    # Following the example from above, new_θ_lims would now be
    #   10-element Vector{Tuple{Int64, Int64}}:
    #   (180, 170)
    #   (170, 150)
    #   (150, 120)
    #   (120, 100)
    #   (100, 90)
    #   (0, 10)
    #   (10, 30)
    #   (30, 60)
    #   (60, 80)
    #   (80, 90)
    angles_to_plot_vert = vcat(eachrow(new_θ_lims_temp)...)
    # @views for i in eachindex(angles_to_plot_vert)
    #     # Find the corresponding indices idx_θ in the simulation grid θ_lims
    #     idx_θ = minimum(angles_to_plot_vert[i]) .<= acosd.(mu_avg(θ_lims)) .<= maximum(angles_to_plot_vert[i])
    #     Ie_plot[i, :, :, :] .= dropdims(sum(Ie[idx_θ, :, :, :]; dims = 1); dims = 1)
    # end
    # for i in eachindex(angles_to_plot_vert)
    #     # Find the corresponding indices idx_θ in the simulation grid θ_lims
    #     idx_θ = minimum(angles_to_plot_vert[i]) .<= acosd.(mu_avg(θ_lims)) .<= maximum(angles_to_plot_vert[i])
    #     @views for j in axes(Ie, 1)[idx_θ]
    #         Ie_plot[i, :, :, :] .+= Ie[j, :, :, :]
    #     end
    # end


    # # Extracts the values, remove the doublets and sort.
    # # Continuing the example from above, this gives us
    # #       [180 170 150 120 100 90 80 60 30 10 0]
    # new_θ_lims_temp = sort(unique(collect(Iterators.flatten(new_θ_lims_temp))); rev=true)

    return Ie_plot
end



#=
This is the function that creates the plot from the data observable. Then when we update the
data, the plot is automatically updated.

It should take Ie as a function of pitch-angle, height and energy as input. But NOT as a
function of time. This because we want to update the time OUTSIDE of the function.
=#
function make_IeztE_3Dzoft_plot(Ie_timeslice::Observable{Array{Float64, 3}},
                                time::Observable{String}, h_atm, E, angles_to_plot, color_limits)

    # Slice the input Ie into its different pitch-angle components
    Ie_streams = Array{Observable}(nothing, length(angles_to_plot))
    for i in eachindex(angles_to_plot)
        Ie_streams[i] = @lift($Ie_timeslice[i, :, :]')
    end

    # Plot (TODO: redo with loop and proper subplots down and up)
    fig = Figure()
    n_row = size(new_θ_lims, 1)
    n_col = size(new_θ_lims, 2)
    for i in axes(new_θ_lims, 1)
        for j in axes(new_θ_lims, 2)
            idx = (i - 1) * n_col  + j # goes along first row, then second row

            ax = Axis(fig[i, j], xscale = log10)
            heatmap!(E, h_atm / 1e3, Ie_streams[idx], colorrange = color_limits, colorscale = log10)

            if i == 1
                ax.title = string(angles_to_plot[i, j][1], " - ", angles_to_plot[i, j][2], "° DOWN")
                ax.xticklabelsvisible = false
            else
                ax.title = string(angles_to_plot[i, j][1], " - ", angles_to_plot[i, j][2], "° UP")
                ax.xlabel = "Energy (eV)"
            end
            if j > 1
                ax.yticklabelsvisible = false
            else
                ax.ylabel = "Height (km)"
            end
        end
    end

    cb = Colorbar(fig[:, end + 1]; limits = color_limits, )
    Label(fig[0, :], time; tellwidth = false, fontsize=18)

    return fig
end



## Call the function

# First load Ie
# filename = "/run/user/1000/gvfs/sftp:host=revontuli.uit.no/mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_536s_longer_fixed-secondaries/IeFlickering-01.mat"
filename = "data/Visions2/IeFlickering-01.mat"
@time data = matread(filename);
@time Ie_raw = data["Ie_ztE"]; # size of [n_mu x nz, nt, nE]

μ_scatterings = data["mu_scatterings"]
μ_lims = data["mu_lims"]
t_run = data["t_run"]
E = data["E"]
dE = diff(E); dE = [dE; dE[end]];
h_atm = data["h_atm"]

@time Ie = restructure_Ie(Ie_raw, μ_lims, h_atm, t_run, E); # size [n_mu, nz, nt, nE]

θ_lims = acosd.(μ_lims)
angles_to_plot = [(0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90);  # DOWN
                 (0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90)]  # UP

@time Ie_plot = restructure_streams_of_Ie(Ie, θ_lims, angles_to_plot);

# Ie_plot_save = copy(Ie_plot);
Ie_plot == Ie_plot_save

Ie_timeslice = Observable(Ie_plot[:, :, 1, :])
time = Observable("$(t_run[1]) s")

global_min, global_max = extrema(Ie_plot)
# global_min, global_max = extrema(Ie)
global_min = 1e0

f = make_IeztE_3Dzoft_plot(Ie_timeslice, time, angles_to_plot, (global_min, global_max))
display(f)














## Test the animation
@time for i_t in 1:n_t
    Ie_timeslice[] .= Ie_plot[:, :, i_t, :]
    time[] = "$(t_run[i_t]) s"
    notify(Ie_timeslice)
    sleep(0.01)
end

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
# Main function
function animate_IeztE_3Dzoft(filename, angles_to_plot)
    # Read the data
    println("Load the data.")

    data = matread(filename);
    Ie_raw = data["Ie_ztE"]; # size of [n_mu x nz, nt, nE]
    # μ_scatterings = data["mu_scatterings"]
    μ_lims = data["mu_lims"]
    t_run = data["t_run"]
    E = data["E"]
    h_atm = data["h_atm"]

    dE = diff(E); dE = [dE; dE[end]]
    θ_lims = acosd.(μ_lims)

    println("Restructure from 3D to 4D array.")
    # Restructure from [n_mu x nz, nt, nE]  to [n_mu, nz, nt, nE]
    Ie = restructure_Ie(Ie_raw, μ_lims, h_atm, t_run, E) # size [n_mu, nz, nt, nE]

    println("Merge the streams to match the angles to plot.")
    # Merge the streams to angles_to_plot and
    Ie_plot = restructure_streams_of_Ie(Ie, θ_lims, angles_to_plot, dE)

    println("Convert from #e-/m²/s to #e-/m²/s/eV/ster.")
    # Convert from #e-/m²/s to #e-/m²/s/eV/ster
    # Restructure `angles_to_plot` from a 2D array to a 1D vector, by concatenating the rows.
    # example: [DOWN1 DOWN2 DOWN3; UP1 UP2 UP3] becomes [DOWN1 DOWN2 DOWN3 UP1 UP2 UP3]
    angles_to_plot_vert = vcat(eachrow(angles_to_plot)...)
    for i in eachindex(angles_to_plot_vert)
        Ie_plot[i] = Ie_plot[i] ./ beam_weight(angles_to_plot_vert[i]) ./ reshape(dE, (1, 1, 1, :))
    end

    global_min, global_max = extrema(Ie_plot)

    println("Create the plot.")
    Ie_timeslice = Observable(Ie_plot[:, :, 1, :])
    time = Observable("$(t_run[1]) s")
    f = make_IeztE_3Dzoft_plot(Ie_timeslice, time, angles_to_plot, (global_min, global_max))

    println("Animate.")
    @time for i_t in 1:n_t
        Ie_timeslice[] .= Ie_plot[:, :, i_t, :]
        time[] = "$(t_run[i_t]) s"
        notify(Ie_timeslice)
        sleep(0.01)
    end

    return nothing
end


##
filename = "data/Visions2/IeFlickering-01.mat"

angles_to_plot = [(0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90);  # DOWN
                 (0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90)]  # UP

animate_IeztE_3Dzoft(filename, angles_to_plot)
