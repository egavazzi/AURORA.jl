#=
Functions to animate Ie(z, t, E).

So far we have
- animate_IeztE_3Dzoft: Ie as a heatmap over height and energy, animation in time.


=#
using CairoMakie
using MAT
using Printf


# Main function
function animate_IeztE_3Dzoft(directory_to_process, angles_to_plot, color_limits; plot_Ietop = false)
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
    Ie = restructure_Ie_from_3D_to_4D(Ie_raw, μ_lims, h_atm, t_run, E) # size [n_mu, nz, nt, nE]
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
    println("Animation will be saved at $video_filename.")
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
                Ie = restructure_Ie_from_3D_to_4D(Ie_raw, μ_lims, h_atm, t_run, E) # size [n_μ, nz, nt, nE]
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
