#=
Functions to animate Ie(z, t, E).

So far we have
- animate_Ie_in_time: Ie as a heatmap over height and energy, animation in time.


=#
using Makie
using NCDatasets: NCDataset
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


# Main function
function AURORA.animate_Ie_in_time(directory_to_process;
                                   angles_to_plot = nothing,
                                   colorrange = nothing,
                                   save_to_file = true,
                                   plot_input = false,
                                   input_angle_cone = [170, 180],
                                   dt_steps = 1,
                                   framerate = 30,
                                   max_bytes = 512 * 1024^2)
    ## Resolve the directory path
    full_path_to_directory = abspath(directory_to_process)
    if !isdir(full_path_to_directory)
        error("Directory not found: '$directory_to_process'")
    end
    if dt_steps < 1
        error("`dt_steps` must be an integer greater than or equal to 1.")
    end

    ## Create file name
    if save_to_file
        video_filename = plot_input ? "animation_with_precipitation.mp4" : "animation.mp4"
        full_video_filename = AURORA.rename_if_exists(joinpath(full_path_to_directory, video_filename))
        println(styled"{bold:The animation will be saved at $video_filename.}")
    end

    ## Load the coordinates
    print("(1/1) [$(Dates.format(now(), "HH:MM:SS"))] Load coordinates... ")
    coord = AURORA.load_coordinates(full_path_to_directory)
    μ_lims    = coord.μ_lims
    t_run     = coord.t
    E_centers = coord.E_centers
    z         = coord.h_atm
    ΔE        = coord.ΔE

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

    # Subsample the time axis
    t_sub = t_run[1:dt_steps:length(t_run)]

    # Compute default colorrange if not provided: a streaming pass over the file to
    # find the global maximum of the processed flux.
    if isnothing(colorrange)
        print("scan colorrange... ")
        max_val = -Inf
        AURORA.foreach_Ie_time_chunk(full_path_to_directory; max_bytes) do Ie_chunk, t_range
            Ie_block, g = process_animation_chunk(Ie_chunk, t_range, dt_steps,
                                                  θ_lims, angles_to_plot, scale_μ, ΔE)
            isempty(g) && return
            max_val = max(max_val, maximum(filter(!isnan, Ie_block)))
        end
        colorrange = (max_val / 1e4, max_val)
    end

    # Load input flux if requested
    input_struct = if plot_input
        nc_path = joinpath(full_path_to_directory, "simulation_data.nc")
        NCDataset(nc_path, "r") do ds
            if !haskey(ds, "Ie_input")
                error("Input flux (`Ie_input`) not found in simulation_data.nc. " *
                      "This file may have been produced by an older version of AURORA.")
            end
            Ietop_raw = Array(ds["Ie_input"])  # [n_μ, n_t_input, n_E]
            t_top_raw = Array(ds["time_input"])
            t_top = [t_top_raw; t_top_raw[end] + diff(t_top_raw)[end]] .- t_top_raw[1]
            idx_θ = vec(input_angle_cone[1] .<= abs.(acosd.(mu_avg(θ_lims))) .<= input_angle_cone[2])
            Ω_beam = beam_weight([input_angle_cone[1], input_angle_cone[2]])
            data_input = dropdims(sum(Ietop_raw[idx_θ, :, :]; dims=1); dims=1) ./
                         Ω_beam ./ ΔE' .* E_centers'
            (; bool = plot_input, t_top, data_input, input_angle_cone)
        end
    else
        (; bool = false, t_top = nothing, data_input = nothing, input_angle_cone = nothing)
    end

    # Create the figure, already containing the first frame
    print("create figure... ")
    first_slice = NCDataset(joinpath(full_path_to_directory, "simulation_data.nc"), "r") do ds
        ds["Ie"].var[:, :, 1:1, :]
    end
    first_frame, = process_animation_chunk(first_slice, 1:1, dt_steps,
                                           θ_lims, angles_to_plot, scale_μ, ΔE)
    Ie_timeslice = Observable(@view first_frame[:, :, 1, :])
    time_obs = Observable(@sprintf("%.3f s", t_sub[1]))
    fig = make_Ie_in_time_plot(Ie_timeslice, time_obs, z, E_centers, angles_to_plot, colorrange, input_struct)
    display(fig)

    # Animate
    if save_to_file
        io = VideoStream(fig; framerate=framerate, visible=true)
    end

    println("animate.")
    # Stream the flux again, recording one frame per subsampled time step. Only one
    # processed time-chunk is held in memory at a time.
    AURORA.foreach_Ie_time_chunk(full_path_to_directory; max_bytes) do Ie_chunk, t_range
        Ie_block, g = process_animation_chunk(Ie_chunk, t_range, dt_steps,
                                              θ_lims, angles_to_plot, scale_μ, ΔE)
        for (k, gk) in enumerate(g)
            Ie_timeslice[] = @view Ie_block[:, :, k, :]
            time_obs[] = @sprintf("%.3f s", t_run[gk])
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


"""
    process_animation_chunk(Ie_chunk, t_range, dt_steps, θ_lims, angles_to_plot, scale_μ, ΔE)

Process one streamed flux chunk for animation: select the subsampled time steps
(`1:dt_steps:end` on the global time axis) that fall inside `t_range`, restructure the beams
into the panels of `angles_to_plot`, and convert to differential flux (#e⁻/m²/s/eV/ster).

Returns `(Ie_block, g)` where `Ie_block` is `[n_z, n_panels, length(g), n_E]` and `g` are the
corresponding global time indices (empty if no subsampled step falls inside `t_range`).
`Ie_block` is a fresh array, safe to keep after `Ie_chunk` (a reused buffer) is overwritten.
"""
function process_animation_chunk(Ie_chunk, t_range, dt_steps, θ_lims, angles_to_plot,
                                 scale_μ, ΔE)
    g = [t for t in t_range if (t - 1) % dt_steps == 0]
    isempty(g) && return zeros(0, 0, 0, 0), g
    Ie_sub = Ie_chunk[:, :, g .- (first(t_range) - 1), :]
    Ie_block = AURORA.restructure_streams_of_Ie(Ie_sub, θ_lims, angles_to_plot)
    n_z, n_panels, n_t, n_E = size(Ie_block)
    @tturbo for i_E in 1:n_E, i_t in 1:n_t, i_μ in 1:n_panels, i_z in 1:n_z
        Ie_block[i_z, i_μ, i_t, i_E] *= scale_μ[i_μ] / ΔE[i_E]
    end
    return Ie_block, g
end
