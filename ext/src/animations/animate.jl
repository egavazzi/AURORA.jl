#=
Functions to animate Ie(z, t, E).

So far we have
- animate_IeztE_3Dzoft: Ie as a heatmap over height and energy, animation in time.


=#
using Makie
using MAT: matread
using Printf: @sprintf
using Term: @bold


# Main function
"""
    animate_IeztE_3Dzoft(directory_to_process, angles_to_plot, color_limits; plot_Ietop = false)

Plot a heatmap of Ie over height and energy, and animate it in time. It will load the
result files one by one. The animation will be saved as a .mp4 file under the
`directory_to_process`.

# Example
```julia-repl
julia> directory_to_process = "Visions2/Alfven_475s";
julia> angles_to_plot = [(0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90);   # DOWN
                         (0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90)];  # UP
julia> color_limits = (1e5, 1e9);
julia> animate_IeztE_3Dzoft(directory_to_process, angles_to_plot, color_limits; plot_Ietop = true)
```

Note that the first row of `angles_to_plot` gives the angles for the down-flux, and the
second row gives the angles for the up-flux. The angles are given from 0° to 90°, and the
function takes care of merging the beams together while keeping the correct units.

*Important*: Right now, the function only supports an `angles_to_plot` with **maximum two rows**.
You can also use only one row which will plot only the down-flux.

*Important*: The limits of the `angles_to_plot` needs to match existing limits of the beams
used in the simulation. E.g. if `θ_lims = 180:-10:0` was used in the simulation, `(30, 60)`
will be fine as 30° and 60° exist as limits, but `(35, 60)` will not as 35° does not exist
as a limit.

# Arguments
- `directory_to_process`: directory to process
- `angles_to_plot`: limits of the angles to plot
- `color_limits`: limits for the colormap/colorbar

# Keyword Arguments
- `plot_Ietop = false`: if true, also plots the precipitating Ie at the top of the
                        ionosphere by loading it from the file `Ie_top.mat`
"""
function AURORA.animate_IeztE_3Dzoft(directory_to_process, angles_to_plot, color_limits; plot_Ietop = false)
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
    θ_lims = round.(acosd.(μ_lims))

    # Restructure from [n_mu x nz, nt, nE]  to [n_mu, nz, nt, nE]
    println("Restructure from 3D to 4D array.")
    Ie = AURORA.restructure_Ie_from_3D_to_4D(Ie_raw, μ_lims, h_atm, t_run, E) # size [n_mu, nz, nt, nE]
    # Merge the streams to angles_to_plot
    println("Merge the streams to match the angles to plot.")
    Ie_plot = AURORA.restructure_streams_of_Ie(Ie, θ_lims, angles_to_plot)
    # Restructure `angles_to_plot` from a 2D array to a 1D vector, by concatenating the rows.
    angles_to_plot_vert = vcat(eachrow(angles_to_plot)...) #
    # Convert from #e-/m²/s to #e-/m²/s/eV/ster
    println("Convert from #e-/m²/s to #e-/m²/s/eV/ster.")
    for i_μ in eachindex(angles_to_plot_vert)
        Ie_plot[i_μ, :, :, :] = Ie_plot[i_μ, :, :, :] ./ beam_weight(angles_to_plot_vert[i_μ]) ./ reshape(dE, (1, 1, :))
    end

    # Load Ietop if requested
    if plot_Ietop
        println("Load Ietop.")
        Ietop_file = find_Ietop_file(full_path_to_directory)
        data = matread(Ietop_file)
        Ietop = data["Ie_total"]
        t_top = data["t_top"]; t_top = [t_top; t_top[end] + diff(t_top)[end]] .- t_top[1]
        angle_cone = [90 180] # angle for the cone of precipitation to plot
        idx_θ = vec(angle_cone[1] .<= abs.(acosd.(mu_avg(θ_lims))) .<= angle_cone[2])
        BeamW = beam_weight([angle_cone[1], angle_cone[2]])
        data_Ietop = dropdims(sum(Ietop[idx_θ, :, :]; dims=1); dims=1) ./ BeamW ./ dE' .* E' # in eV/m²/s/eV/ster
        Ietop_struct = (bool = plot_Ietop, t_top = t_top, data_Ietop = data_Ietop)
        video_filename = rename_if_exists(joinpath(full_path_to_directory, "animation_with_precipitation.mp4"))
    else
        Ietop_struct = (bool = false, t_top = nothing, data_Ietop = nothing)
        video_filename = rename_if_exists(joinpath(full_path_to_directory, "animation.mp4"))
    end

    # Create the figure
    println("Create the plot.")
    Ie_timeslice = Observable(Ie_plot[:, :, 1, :])
    time = Observable("$(t_run[1]) s")
    fig = make_IeztE_3Dzoft_plot(Ie_timeslice, time, h_atm, E, angles_to_plot, color_limits, Ietop_struct)
    display(fig)

    # Animate
    println(@bold "The animation will be saved at $video_filename.")
    println("Animate.")
    n_t = length(t_run) # number of timesteps per file
    record(fig, video_filename) do io
        for i_file in 1:n_files
            # First file (already loaded)
            if i_file == 1
                i_t = 1
                while i_t <= n_t
                    # string(Makie.current_backend()) == "GLMakie" || isopen(fig.scene) || error("Figure is closed") # exit if GLMakie window is closed
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
                print("($i_file/$n_files) Load data... ")
                data = matread(files_to_process[i_file]);
                Ie_raw = data["Ie_ztE"]; # size of [n_μ x nz, nt, nE]
                t_run = data["t_run"]
                # Restructure from [n_mu x nz, nt, nE]  to [n_mu, nz, nt, nE]
                Ie = AURORA.restructure_Ie_from_3D_to_4D(Ie_raw, μ_lims, h_atm, t_run, E) # size [n_μ, nz, nt, nE]
                # Merge the streams to angles_to_plot
                Ie_plot = AURORA.restructure_streams_of_Ie(Ie, θ_lims, angles_to_plot) # size [n_μ_new, nz, nt, nE]
                # Convert from #e-/m²/s to #e-/m²/s/eV/ster
                for i_μ in eachindex(angles_to_plot_vert)
                    Ie_plot[i_μ, :, :, :] = Ie_plot[i_μ, :, :, :] ./ beam_weight(angles_to_plot_vert[i_μ]) ./ reshape(dE, (1, 1, :))
                end
                # Animate
                println("Animate.")
                i_t = 2 # we skip the first timestamp as it is the same as the last element from the previous file.
                while i_t <= n_t
                    # string(Makie.current_backend()) == "GLMakie" || isopen(fig.scene) || error("Figure is closed") # exit if GLMakie window is closed
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
