#=
Functions to animate Ie(z, t, E).

So far we have
- animate_Ie_in_time: Ie as a heatmap over height and energy, animation in time.


=#
using Makie
using MAT: matread
using Printf: @sprintf
using Term: @bold
using Dates


"""
    generate_default_angles(θ_lims)

Generate a default `angles_to_plot` matrix from the simulation's θ_lims.
Returns a 2-row matrix where row 1 is down-flux (angles >= 90°) and row 2 is up-flux
(angles <= 90°). Beams that cross 90° are included in the down-flux row.
Angles are in the range 0-180° where 180° is field-aligned down and 0° is field-aligned up.
"""
function generate_default_angles(θ_lims)
    # Sort θ_lims from 180 to 0
    θ_sorted = sort(collect(θ_lims); rev=true)

    # Build all beams as consecutive pairs
    all_beams = [(θ_sorted[i], θ_sorted[i+1]) for i in 1:(length(θ_sorted)-1)]

    # Separate into down (includes beams that touch or cross 90°) and up (fully below 90°)
    angles_down = Union{Tuple{Float64, Float64}, Nothing}[b for b in all_beams if b[1] >= 90]  # beam starts at or above 90°
    angles_up = Union{Tuple{Float64, Float64}, Nothing}[b for b in all_beams if b[1] < 90]     # beam starts below 90°

    # Pad shorter row with nothing to match dimensions
    n_cols = max(length(angles_down), length(angles_up))
    while length(angles_down) < n_cols
        push!(angles_down, nothing)
    end
    while length(angles_up) < n_cols
        push!(angles_up, nothing)
    end

    # Build matrix: row 1 = DOWN, row 2 = UP
    return vcat(permutedims(angles_down), permutedims(angles_up))
end


"""
    compute_default_colorrange(Ie_plot)

Compute a default colorrange based on the maximum value of the data.
Returns (max_value / 1e4, max_value) to span 4 orders of magnitude.
"""
function compute_default_colorrange(Ie_plot)
    max_val = maximum(filter(!isnan, Ie_plot))
    return (max_val / 1e4, max_val)
end


# Main function
function AURORA.animate_Ie_in_time(directory_to_process;
                                     angles_to_plot = nothing,
                                     colorrange = nothing,
                                     save_to_file = true,
                                     plot_Ietop = false,
                                     Ietop_angle_cone = [170, 180])
    ## Resolve the directory path
    full_path_to_directory = abspath(directory_to_process)
    if !isdir(full_path_to_directory)
        error("Directory not found: '$directory_to_process'")
    end

    ## Find the files to process
    files = readdir(full_path_to_directory, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]
    n_files = length(files_to_process)

    ## Create file name
    if save_to_file
        video_filename = plot_Ietop ? "animation_with_precipitation.mp4" : "animation.mp4"
        full_video_filename = rename_if_exists(joinpath(full_path_to_directory, video_filename))
        println(@bold "The animation will be saved at $video_filename.")
    end

    ## Load data from the first file
    print("(1/$n_files) [$(Dates.format(now(), "HH:MM:SS"))] Load data... ")
    data = matread(files_to_process[1])
    Ie_raw = data["Ie_ztE"]  # size of [n_mu x nz, nt, nE]
    μ_lims = vec(data["mu_lims"])
    t_run = data["t_run"]
    E = data["E"]
    h_atm = data["h_atm"]

    dE = diff(E); dE = [dE; dE[end]]
    θ_lims = round.(acosd.(μ_lims))

    # Use default angles_to_plot if not provided
    if isnothing(angles_to_plot)
        angles_to_plot = generate_default_angles(θ_lims)
    end

    # Restructure `angles_to_plot` from a 2D array to a 1D vector
    angles_to_plot_vert = vec(permutedims(angles_to_plot))

    # Helper function to process Ie data
    function process_Ie(Ie_raw, t_run)
        Ie = AURORA.restructure_Ie_from_3D_to_4D(Ie_raw, μ_lims, h_atm, t_run, E)
        Ie_plot = AURORA.restructure_streams_of_Ie(Ie, θ_lims, angles_to_plot)
        # Convert from #e-/m²/s to #e-/m²/s/eV/ster
        for i_μ in eachindex(angles_to_plot_vert)
            isnothing(angles_to_plot_vert[i_μ]) && continue  # skip empty panels
            θ1, θ2 = angles_to_plot_vert[i_μ]
            Ie_plot[i_μ, :, :, :] ./= beam_weight([θ1, θ2]) .* reshape(dE, (1, 1, :))
        end
        return Ie_plot
    end

    Ie_plot = process_Ie(Ie_raw, t_run)

    # Compute default colorrange if not provided
    if isnothing(colorrange)
        colorrange = compute_default_colorrange(Ie_plot)
    end

    # Load Ietop if requested
    if plot_Ietop
        Ietop_file = find_Ietop_file(full_path_to_directory)
        data = matread(Ietop_file)
        Ietop = data["Ie_total"]
        t_top = data["t_top"]; t_top = [t_top; t_top[end] + diff(t_top)[end]] .- t_top[1]
        idx_θ = vec(Ietop_angle_cone[1] .<= abs.(acosd.(mu_avg(θ_lims))) .<= Ietop_angle_cone[2])
        BeamW = beam_weight([Ietop_angle_cone[1], Ietop_angle_cone[2]])
        data_Ietop = dropdims(sum(Ietop[idx_θ, :, :]; dims=1); dims=1) ./ BeamW ./ dE' .* E' # in eV/m²/s/eV/ster
        Ietop_struct = (; bool = plot_Ietop, t_top, data_Ietop, Ietop_angle_cone)
    else
        Ietop_struct = (; bool = false, t_top = nothing, data_Ietop = nothing, Ietop_angle_cone = nothing)
    end

    # Create the figure
    print("create figure... ")
    Ie_timeslice = Observable(Ie_plot[:, :, 1, :])
    time = Observable("$(t_run[1]) s")
    fig = make_Ie_in_time_plot(Ie_timeslice, time, h_atm, E, angles_to_plot, colorrange, Ietop_struct)
    display(fig)

    # Animate
    framerate = 30
    if save_to_file
        io = VideoStream(fig; framerate=framerate, visible=true)
    end

    for i_file in 1:n_files
        # Load data (first file already loaded)
        if i_file > 1
            print("($i_file/$n_files) [$(Dates.format(now(), "HH:MM:SS"))] Load data... ")
            data = matread(files_to_process[i_file])
            Ie_raw = data["Ie_ztE"]
            t_run = data["t_run"]
            Ie_plot = process_Ie(Ie_raw, t_run)
        end

        # Animate (skip first timestep for files after the first, as it duplicates the
        # previous file's last frame)
        i_t_start = i_file == 1 ? 1 : 2
        println("animate.")
        for i_t in i_t_start:length(t_run)
            Ie_timeslice[] .= Ie_plot[:, :, i_t, :]
            time[] = @sprintf("%.3f s", t_run[i_t])
            notify(Ie_timeslice)
            if save_to_file
                recordframe!(io)
            end
            sleep(0.005)
        end
    end

    if save_to_file
        print("\nSaving animation... ")
        save(full_video_filename, io)
        println("done!")
    end

    return fig
end
