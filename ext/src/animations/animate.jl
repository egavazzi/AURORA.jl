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
                                   plot_input = false,
                                   input_angle_cone = [170, 180],
                                   dt_steps = 1,
                                   framerate = 30)
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

    ## Load simulation data
    print("(1/1) [$(Dates.format(now(), "HH:MM:SS"))] Load data... ")
    result = AURORA.read_simulation_nc(full_path_to_directory)
    Ie       = result.Ie          # [n_z, n_μ, n_t, n_E]
    μ_lims   = result.mu_lims
    t_run    = result.t
    E_centers = result.E_centers
    z        = result.h_atm
    ΔE       = result.dE

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
    t_indices = 1:dt_steps:length(t_run)
    t_sub     = t_run[t_indices]
    Ie_sub    = Ie[:, :, t_indices, :]

    # Restructure and normalise beams
    Ie_plot = AURORA.restructure_streams_of_Ie(Ie_sub, θ_lims, angles_to_plot)
    # Ie_plot is [n_z, n_panels, n_t, n_E] — convert to #e-/m²/s/eV/ster
    n_z, n_panels, n_t_sub, n_E = size(Ie_plot)
    @tturbo for i_E in 1:n_E, i_t in 1:n_t_sub, i_μ in 1:n_panels, i_z in 1:n_z
        Ie_plot[i_z, i_μ, i_t, i_E] *= scale_μ[i_μ] / ΔE[i_E]
    end

    # Compute default colorrange if not provided
    if isnothing(colorrange)
        colorrange = compute_default_colorrange(Ie_plot)
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

    # Create the figure
    print("create figure... ")
    Ie_timeslice = Observable(@view Ie_plot[:, :, 1, :])
    time_obs = Observable("$(t_sub[1]) s")
    fig = make_Ie_in_time_plot(Ie_timeslice, time_obs, z, E_centers, angles_to_plot, colorrange, input_struct)
    display(fig)

    # Animate
    if save_to_file
        io = VideoStream(fig; framerate=framerate, visible=true)
    end

    println("animate.")
    for i_t in eachindex(t_sub)
        Ie_timeslice[] = @view Ie_plot[:, :, i_t, :]
        time_obs[] = @sprintf("%.3f s", t_sub[i_t])
        if save_to_file
            recordframe!(io)
        end
        sleep(0.01)
    end

    if save_to_file
        print("\nSaving animation... ")
        save(full_video_filename, io)
        println("done!")
    end

    return fig
end
