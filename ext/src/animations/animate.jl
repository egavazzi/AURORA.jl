#=
Functions to animate Ie(z, t, E).

So far we have
- animate_Ie_in_time: Ie as a heatmap over height and energy, animation in time.


=#
using Makie
using MAT: matread
using Printf: @sprintf
using StyledStrings: @styled_str
using Dates
using LoopVectorization: @tturbo


"""
    generate_default_angles(θ_lims)

Generate a default `angles_to_plot` matrix from the simulation's θ_lims.
Returns a 2-row matrix where row 1 is down-flux (negative average μ) and row 2 is
up-flux (positive average μ).
Angles are in the range 0-180° where 180° is field-aligned down and 0° is field-aligned up.
"""
function generate_default_angles(θ_lims)
    # Sort θ_lims from 180 to 0
    θ_sorted = sort(collect(θ_lims); rev=true)

    # Build all beams as consecutive pairs
    all_beams = [(θ_sorted[i], θ_sorted[i+1]) for i in 1:(length(θ_sorted)-1)]
    all_beams_μ = mu_avg(θ_sorted)

    # Separate into down and up using the sign of the average pitch-angle cosine.
    # Row 1: most field-aligned down-going (μ ≈ -1) to near-perpendicular (μ ≈ 0).
    # Row 2: most field-aligned up-going (μ ≈ 1) to near-perpendicular (μ ≈ 0).
    idx_down = sort([i for i in eachindex(all_beams_μ) if all_beams_μ[i] < 0]; by = i -> all_beams_μ[i])
    idx_up = sort([i for i in eachindex(all_beams_μ) if all_beams_μ[i] > 0]; by = i -> all_beams_μ[i], rev = true)

    angles_down = Union{Tuple{Float64, Float64}, Nothing}[all_beams[i] for i in idx_down]
    angles_up = Union{Tuple{Float64, Float64}, Nothing}[(all_beams[i][2], all_beams[i][1]) for i in idx_up]

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
                                     Ietop_angle_cone = [170, 180],
                                     dt_steps = 1)
    ## Resolve the directory path
    full_path_to_directory = abspath(directory_to_process)
    if !isdir(full_path_to_directory)
        error("Directory not found: '$directory_to_process'")
    end
    if dt_steps < 1
        error("`dt_steps` must be an integer greater than or equal to 1.")
    end

    ## Find the files to process
    files = readdir(full_path_to_directory, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]
    n_files = length(files_to_process)

    ## Create file name
    if save_to_file
        video_filename = plot_Ietop ? "animation_with_precipitation.mp4" : "animation.mp4"
        full_video_filename = rename_if_exists(joinpath(full_path_to_directory, video_filename))
        println(styled"{bold:The animation will be saved at $video_filename.}")
    end

    ## Load data from the first file
    print("(1/$n_files) [$(Dates.format(now(), "HH:MM:SS"))] Load data... ")
    data = matread(files_to_process[1])
    Ie_raw = data["Ie_ztE"]  # size of [n_mu x nz, nt, nE]
    μ_lims = vec(data["mu_lims"])
    t_run = data["t_run"]
    E_centers = data["E_centers"]
    z = data["h_atm"]

    ΔE = data["dE"]
    θ_lims = round.(acosd.(μ_lims))

    # Use default angles_to_plot if not provided
    if isnothing(angles_to_plot)
        angles_to_plot = generate_default_angles(θ_lims)
    end

    # Restructure `angles_to_plot` from a 2D array to a 1D vector
    angles_to_plot_vert = vec(permutedims(angles_to_plot))

    # Precompute solid-angle normalization for non-empty panels.
    scale_μ = zeros(Float64, length(angles_to_plot_vert))
    for i_μ in eachindex(angles_to_plot_vert)
        θ_bin = angles_to_plot_vert[i_μ]
        isnothing(θ_bin) && continue
        θ1, θ2 = θ_bin
        scale_μ[i_μ] = 1 ./ beam_weight([θ1, θ2])[1]
    end

    Ie_plot_buffer = Ref{Union{Nothing, Array{Float64, 4}}}(nothing)

    # Helper function to process Ie data
    function process_Ie(Ie_raw, t_run)
        Ie = AURORA.restructure_Ie_from_3D_to_4D(Ie_raw, μ_lims, z, t_run, E_centers)
        buffer_size = (size(Ie, 1), length(angles_to_plot), size(Ie, 3), size(Ie, 4))
        if isnothing(Ie_plot_buffer[]) || size(Ie_plot_buffer[]) != buffer_size
            Ie_plot_buffer[] = Array{Float64, 4}(undef, buffer_size...)
        end

        Ie_plot = AURORA.restructure_streams_of_Ie!(Ie_plot_buffer[], Ie, θ_lims, angles_to_plot)
        # Convert from #e-/m²/s to #e-/m²/s/eV/ster
        for i_μ in eachindex(angles_to_plot_vert)
            isnothing(angles_to_plot_vert[i_μ]) && continue  # skip empty panels
            @tturbo for i_E in eachindex(E_centers), i_t in eachindex(t_run), i_z in eachindex(z)
                Ie_plot[i_z, i_μ, i_t, i_E] *= scale_μ[i_μ] / ΔE[i_E]
            end
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
        Ω_beam = beam_weight([Ietop_angle_cone[1], Ietop_angle_cone[2]])
        data_Ietop = dropdims(sum(Ietop[idx_θ, :, :]; dims=1); dims=1) ./ Ω_beam ./ ΔE' .* E_centers' # in eV/m²/s/eV/ster
        Ietop_struct = (; bool = plot_Ietop, t_top, data_Ietop, Ietop_angle_cone)
    else
        Ietop_struct = (; bool = false, t_top = nothing, data_Ietop = nothing, Ietop_angle_cone = nothing)
    end

    # Create the figure
    print("create figure... ")
    Ie_timeslice = Observable(@view Ie_plot[:, :, 1, :])
    time = Observable("$(t_run[1]) s")
    fig = make_Ie_in_time_plot(Ie_timeslice, time, z, E_centers, angles_to_plot, colorrange, Ietop_struct)
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
        for i_t in i_t_start:dt_steps:length(t_run)
            Ie_timeslice[] = @view Ie_plot[:, :, i_t, :]
            time[] = @sprintf("%.3f s", t_run[i_t])
            if save_to_file
                recordframe!(io)
            end
            sleep(0.01)
        end
    end

    if save_to_file
        print("\nSaving animation... ")
        save(full_video_filename, io)
        println("done!")
    end

    return fig
end
