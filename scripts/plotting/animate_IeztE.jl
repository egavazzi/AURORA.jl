#=
Functions to animate Ie(z, t, E).

So far we have
- animate_IeztE_3Dzoft: Ie as a heatmap over height and energy, animation in time.


=#
using GLMakie
using CairoMakie
using MAT
using AURORA
using Printf
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
    restructure_streams_of_Ie(Ie, θ_lims, new_θ_lims)

Function that merges the streams of `Ie` that are given over `θ_lims` to fit the
`new_θ_lims` of interest. It can be useful when wanting to merge some streams for plotting.

For example, if we have
```
    θ_lims = [180 160 140 120 100 90 80 60 40 20 0] # simulation
```
and we want to plot with (this should be an array of tuples, but to simplify the comparison we write it as a vector here)
```
    new_θ_lims = [180 160 120 100 90 80 60 40 20 0] # to plot
```
the function will merge the streams (160°-140°) and (140°-120°) together into a new stream with limits (160°-120°).

# Calling
`Ie_plot = restructure_streams_of_Ie(Ie, θ_lims, new_θ_lims)`

# Arguments
- `Ie`: array of electron flux with pitch-angle limits `θ_lims`. Of shape [n\\_μ, n\\_z, n\\_t, n\\_E].
- `θ_lims`: pitch-angle limits. Usually a vector or range.
- `new_θ_lims`: new pitch-angle limits. Given as an array of tuples with two rows, for example:
```
            julia> new_θ_lims = [(0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90);  # DOWN
                                 (0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90)]  # UP
```

# Returns
- `Ie_plot`: array of electron flux with the new pitch-angle limits `new_θ_lims`. Of shape
             [n\\_μ\\_new, n\\_z, n\\_t, n\\_E], where n\\_μ\\_new is the number of streams
             in `new_θ_lims`. The first dimension of `Ie_plot` is sorted such that the
             indices go along the first row of `new_θ_lims`, and then the second row.
             In our example with `new_θ_lims` from above, that would be ``[1 2 3 4 5; 6 7 8 9 10]``.

"""
function restructure_streams_of_Ie(Ie, θ_lims, new_θ_lims)
    # Initialize the new Ie_plot that will contain the restructured streams
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

    # Restructure `new_θ_lims_temp` from a 2D array to a 1D vector, by concatenating the rows.
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
    new_θ_lims_temp = vcat(eachrow(new_θ_lims_temp)...)

    # Restructure to [n_μ_new, n_z, n_t, n_E]
    # Loop over the new_θ_lims streams
    @views for i in eachindex(new_θ_lims_temp)
        # Find the indices of the streams from the simulation that should be merged in the stream new_θ_lims[i].
        idx_θ = axes(Ie, 1)[minimum(new_θ_lims_temp[i]) .<= acosd.(mu_avg(θ_lims)) .<= maximum(new_θ_lims_temp[i])]
        # Loop over these streams and add them into the right stream of Ie_plot.
        for j in idx_θ
            Ie_plot[i, :, :, :] .+= Ie[j, :, :, :]
        end
    end

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
                                time::Observable{String}, h_atm, E, angles_to_plot, color_limits,
                                t_top, data_Ietop)

    # Slice the input Ie into its different pitch-angle components
    Ie_streams = Array{Observable}(nothing, length(angles_to_plot))
    for i in eachindex(angles_to_plot)
        Ie_streams[i] = @lift($Ie_timeslice[i, :, :]')
    end

    # Plot (TODO: redo with loop and proper subplots down and up)
    fig = Figure(size = (1500, 1000), fontsize=20)
    n_row = size(angles_to_plot, 1)
    n_col = size(angles_to_plot, 2)
    ga = fig[1:4, 1:n_col] = GridLayout()
    for i in axes(angles_to_plot, 1)
        for j in axes(angles_to_plot, 2)
            idx = (i - 1) * n_col  + j # goes along first row, then second row

            ax = Axis(ga[i, j], xscale = log10)
            heatmap!(E, h_atm / 1e3, Ie_streams[idx], colorrange = color_limits, colorscale = log10, colormap = :inferno)

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
    Colorbar(fig[:, end + 1]; limits = color_limits, scale = log10, label = "Ie (#e⁻/m²/s/eV/ster)", colormap = :inferno)
    Label(fig[0, 4], time; tellwidth = false, tellheight = false, fontsize=20)

    # Plot Ie precipitating at the top (#TODO: redo this properly)
    gb = fig[0, 1:2] = GridLayout()
    ax_Ietop = Axis(gb[1, 1], yscale = log10, ylabel = " Energy (eV)", xlabel = "t (s)",
                    title = "Incoming flux at the top",
                    # title = string(180 - angle_cone[2], " - ", 180 - angle_cone[1], "° DOWN"),
                    yminorticksvisible = true, yminorticks = IntervalsBetween(9),
                    xticklabelsvisible = true, xminorticksvisible = true,
                    xticksmirrored = true, yticksmirrored = true,
                    limits = ((0, 1), nothing))
    hm = heatmap!(t_top, E, data_Ietop; colormap = :inferno, colorscale=log10, colorrange=(1e6, maximum(data_Ietop)))
    time_float64 = @lift(parse(Float64, $time[1:end-1]))
    vlines!(time_float64, linewidth = 3)
    Colorbar(gb[1, 2], hm; label = "IeE (eV/m²/s/eV/ster)")

    return fig
end




# Main function
function animate_IeztE_3Dzoft(directory_to_process, angles_to_plot, color_limits)
    ## Find the files to process
    full_path_to_directory = pkgdir(AURORA, "data", directory_to_process)
    files = readdir(full_path_to_directory, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]
    n_files = length(files_to_process)

    ## Load data from the first file
    println("Load the data.")
    data = matread(files_to_process[1]);
    Ie_raw = data["Ie_ztE"]; # size of [n_mu x nz, nt, nE]
    μ_lims = vec(data["mu_lims"])
    t_run = data["t_run"]
    E = data["E"]
    h_atm = data["h_atm"]

    dE = diff(E); dE = [dE; dE[end]]
    θ_lims = acosd.(μ_lims)

    # Restructure from [n_mu x nz, nt, nE]  to [n_mu, nz, nt, nE]
    println("Restructure from 3D to 4D array.")
    Ie = restructure_Ie(Ie_raw, μ_lims, h_atm, t_run, E) # size [n_mu, nz, nt, nE]
    # Merge the streams to angles_to_plot and
    println("Merge the streams to match the angles to plot.")
    Ie_plot = restructure_streams_of_Ie(Ie, θ_lims, angles_to_plot)
    # Convert from #e-/m²/s to #e-/m²/s/eV/ster
    println("Convert from #e-/m²/s to #e-/m²/s/eV/ster.")
    # Restructure `angles_to_plot` from a 2D array to a 1D vector, by concatenating the rows.
    angles_to_plot_vert = vcat(eachrow(angles_to_plot)...)
    for i_μ in eachindex(angles_to_plot_vert)
        Ie_plot[i_μ, :, :, :] = Ie_plot[i_μ, :, :, :] ./ beam_weight(angles_to_plot_vert[i_μ]) ./ reshape(dE, (1, 1, :))
    end

    # Create the figure
    println("Create the plot.")
    Ie_timeslice = Observable(Ie_plot[:, :, 1, :])
    time = Observable("$(t_run[1]) s")
    fig = make_IeztE_3Dzoft_plot(Ie_timeslice, time, h_atm, E, angles_to_plot, color_limits, t_top, data_Ietop)
    display(fig)

    # Animate
    println("Animate.")
    n_t = length(t_run) # number of timesteps per file
    record(fig, joinpath(full_path_to_directory, "animation3.mp4")) do io
        for i_file in 1:n_files
            # First file (already loaded)
            if i_file == 1
                i_t = 1
                while i_t <= n_t
                    Makie.current_backend() == GLMakie || isopen(fig.scene) || error("Figure is closed") # exit if GLMakie window is closed
                    Ie_timeslice[] .= Ie_plot[:, :, i_t, :]
                    time[] = @sprintf("%.3f s", t_run[i_t])
                    notify(Ie_timeslice)
                    sleep(0.01)
                    i_t += 1
                    recordframe!(io)
                end
            # Next files
            else
                # Load the data
                println("Load data... $i_file/$n_files")
                data = matread(files_to_process[i_file]);
                Ie_raw = data["Ie_ztE"]; # size of [n_μ x nz, nt, nE]
                t_run = data["t_run"]
                # Restructure from [n_mu x nz, nt, nE]  to [n_mu, nz, nt, nE]
                Ie = restructure_Ie(Ie_raw, μ_lims, h_atm, t_run, E) # size [n_μ, nz, nt, nE]
                # Merge the streams to angles_to_plot
                Ie_plot = restructure_streams_of_Ie(Ie, θ_lims, angles_to_plot) # size [n_μ_new, nz, nt, nE]
                # Convert from #e-/m²/s to #e-/m²/s/eV/ster
                for i_μ in eachindex(angles_to_plot_vert)
                    Ie_plot[i_μ, :, :, :] = Ie_plot[i_μ, :, :, :] ./ beam_weight(angles_to_plot_vert[i_μ]) ./ reshape(dE, (1, 1, :))
                end
                # Animate
                println("Animate.")
                i_t = 2 # we skip the first timestamp as it is the same as the last element from the previous file.
                while i_t <= n_t
                    Makie.current_backend() == GLMakie || isopen(fig.scene) || error("Figure is closed") # exit if GLMakie window is closed
                    Ie_timeslice[] .= Ie_plot[:, :, i_t, :]
                    time[] = @sprintf("%.3f s", t_run[i_t])
                    notify(Ie_timeslice)
                    sleep(0.01)
                    i_t += 1
                    recordframe!(io)
                end
            end
        end
    end

    return nothing
end


## Calling the animate function

directory_to_process = "Visions2/Alfven_475s_longer_fixed-secondaries"
angles_to_plot = [(0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90);  # DOWN
                  (0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90)]  # UP
color_limits = (1e5, 1e9)
animate_IeztE_3Dzoft(directory_to_process, angles_to_plot, color_limits)






## Testing

# First load Ie
# filename = "/run/user/1000/gvfs/sftp:host=revontuli.uit.no/mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_536s_longer_fixed-secondaries/IeFlickering-01.mat"
filename = "data/Visions2/Alfven_475s_longer_fixed-secondaries/IeFlickering-01.mat"

data = matread(filename);
Ie_raw = data["Ie_ztE"]; # size of [n_mu x nz, nt, nE]
μ_scatterings = data["mu_scatterings"]
μ_lims = data["mu_lims"]
t_run = data["t_run"]
E = data["E"]
h_atm = data["h_atm"]
dE = diff(E); dE = [dE; dE[end]];
θ_lims = acosd.(μ_lims)

angles_to_plot = [(0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90);  # DOWN
                 (0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90)]  # UP
Ie = restructure_Ie(Ie_raw, μ_lims, h_atm, t_run, E); # size [n_mu, nz, nt, nE]
Ie_plot = restructure_streams_of_Ie(Ie, θ_lims, angles_to_plot);

Ie_timeslice = Observable(Ie_plot[:, :, 1, :])
time = Observable("$(t_run[1]) s")
f = make_IeztE_3Dzoft_plot(Ie_timeslice, time, h_atm, E, angles_to_plot, (1e5, 1e9), t_top, data_Ietop)
display(f)


## Test the animation
@time for i_t in 1:length(t_run)
    Ie_timeslice[] .= Ie_plot[:, :, i_t, :]
    time[] = "$(t_run[i_t]) s"
    notify(Ie_timeslice)
    sleep(0.01)
end

## Testing recording the animation to a mp4
# with glmakie(~12s for n_t=51)
@time record(f, "time_animation1.mp4", 1:10; framerate=10, backend=GLMakie) do i_t
    Ie_timeslice[] .= Ie_plot[:, :, i_t, :]
    time[] = "$(t_run[i_t]) s"
    notify(Ie_timeslice)
end
# with cairomakie (~110s for n_t=51)
@time record(f, "time_animation2.mp4", 1:10; framerate=10, backend=CairoMakie) do i_t
    Ie_timeslice[] .= Ie_plot[:, :, i_t, :]
    time[] = "$(t_run[i_t]) s"
    notify(Ie_timeslice)
end






## TODO: redo this properly by integrating this in the animate function (and maybe with a bool to choose to plot it or not?)
full_path_to_directory = "data/Visions2/Alfven_475s_longer_fixed-secondaries"
incoming_files = filter(file -> startswith(file, "Ie_incoming_"), readdir(full_path_to_directory))
if length(incoming_files) > 1
    error("More than one file contains incoming flux. This is not normal")
else
    global Ietop_file = joinpath(full_path_to_directory, incoming_files[1])
end

# Read the Ie_top flux
data = matread(Ietop_file)
Ietop = data["Ie_total"]
t_top = data["t_top"]; t_top = [t_top; t_top[end] + diff(t_top)[end]] .- t_top[1]

angle_cone = [90 180] # angle for the cone of precipitation to plot
idx_θ = vec(angle_cone[1] .<= abs.(acosd.(mu_avg(θ_lims))) .<= angle_cone[2])
BeamW = beam_weight([angle_cone[1], angle_cone[2]])
data_Ietop = dropdims(sum(Ietop[idx_θ, :, :]; dims=1); dims=1) ./ BeamW ./ dE' .* E'
