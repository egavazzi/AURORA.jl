using DataInterpolations: PCHIPInterpolation, ExtrapolationType
using LoopVectorization: @tturbo


"""
    v_of_E(E)

Calculate the velocity (in **m/s**) of an electron with energy `E` (in **eV**).

# Calling
`v = v_of_E(E)`

# Input
- `E` : energy in **eV**, can be a scalar, vector, range, ...

# Output
- `v` : velocity in **m/s**
"""
function v_of_E(E)
	v = (2 * qₑ * abs(E) / mₑ) .^ (1/2) .* sign(E);
	return v
end


using HCubature: hcubature
"""
    mu_avg(θ_lims)

Calculate the cosinus of the center of each pitch-angle beams delimited by θ_lims. This is
for isotropically distributed fluxes within each beam, i.e the fluxes are weighted by sin(θ)

# Calling
`μ_center = mu_avg(θ_lims) `

# Inputs
- `θ_lims` : pitch-angle limits *in degrees* of all the beams, range or vector [n_beams + 1]

# Outputs
- `μ_center` : cosine of the center of all the pitch-angle beams, vector [n_beams + 1]
"""
function mu_avg(θ_lims)
    μ_center = zeros(length(θ_lims) - 1)
    for i in eachindex(μ_center)
        μ_center[i] = hcubature(x -> cosd.(x) .* sind.(x), [θ_lims[i]], [θ_lims[i + 1]])[1][1] ./
                      hcubature(x -> sind.(x), [θ_lims[i]], [θ_lims[i + 1]])[1][1]
    end
    return μ_center
end


using QuadGK: quadgk
"""
    beam_weight(θ_lims)

Return the beam weights of the pitch-angle beams delimited by θ_lims.

# Calling
`BeamW = beam_weight(θ_lims) `

# Inputs
- `θ_lims` : pitch-angle limits of all the beams, range or vector [n_beams + 1]

# Outputs
- `BeamW` : solid angle of each pitch-angle beams (ster), vector [n_beams]
"""
function beam_weight(θ_lims)
    BeamW = Vector{Float64}(undef, length(θ_lims) - 1)
    for i_μ in eachindex(BeamW)
        BeamW[i_μ] = 2 * π * abs(quadgk(sin, deg2rad(θ_lims[i_μ]), deg2rad(θ_lims[i_μ + 1]))[1])
    end
    return BeamW
end


"""
    rename_if_exists(savefile)

This function takes a string as an input. If a file or folder with that name *does not* exist,
it returns the same string back. But if the folder or file already exists, it appends a
number between parenthesis to the name string.

For example, if the folder `foo/` already exist and `"foo"` is given as input, the
function will return a string `"foo(1)"` as an output. Similarly, if a file `foo.txt`
already exists and `"foo.txt"` is given as input, the function will return a string
`"foo(1).txt"`. If the file `foo(1).txt` also already exist, the function will return a string
`"foo(2).txt"`, etc...

The function should support all types of extensions.

# Calling
`newsavefile = rename_if_exists(savefile)`
"""
function rename_if_exists(savefile)
    # Check if path exists (either as file or directory)
    if !ispath(savefile)
        return savefile
    end

    # Remove trailing slash (if any)
    savefile_clean = rstrip(savefile, '/')

    # Split filename and extension
    dir, name = splitdir(savefile_clean)
    name_without_ext, ext = splitext(name)

    # Find next available counter
    counter = 1
    while true
        new_name = "$(name_without_ext)($(counter))$(ext)"
        new_path = joinpath(dir, new_name)
        ispath(new_path) || return new_path
        counter += 1
    end
end


"""
    smooth_transition(x, x_start = 0.0, x_end = 1.0)

Create a smooth transition from 0 to 1 over the interval `[x_start, x_end]` using a
C∞-smooth (infinitely differentiable) function.

- Returns 0 for `x < x_start`
- Returns 1 for `x > x_end`
- Smoothly transitions from 0 to 1 in between, with all derivatives continuous

# Arguments
- `x`: The input value at which to evaluate the transition
- `x_start`: Start of the transition interval (default: 0.0)
- `x_end`: End of the transition interval (default: 1.0)

# Returns
- A value between 0 and 1 representing the smooth transition
"""
function smooth_transition(x, x_start = 0.0, x_end = 1.0)
    # If the interval has zero width, return step function
    if iszero(x_end - x_start)
        return x < x_start ? 0.0 : 1.0
    end

    # Normalize x to the interval [0, 1]
    x_normalized = (x - x_start) / (x_end - x_start)

    # Smoothing function: exp(-1/x) for x > 0, otherwise 0
    Ψ(x) = x <= 0 ? 0.0 : exp(-1 / x)
    smoothing_lower = Ψ(x_normalized)      # Activates for x > x_start
    smoothing_upper = Ψ(1 - x_normalized)  # Activates for x < x_end

    # Combine the parts to create smooth transition
    transition_value = smoothing_lower / (smoothing_lower + smoothing_upper)

    return transition_value
end


# Function to restructure the matrix from 3D [n_mu x nz, nt, nE] to 4D [nz, n_mu, nt, nE]
function restructure_Ie_from_3D_to_4D(Ie_raw, μ_lims, z, t_run, E_centers)
    n_μ = length(μ_lims) - 1
    n_z = length(z)
    n_t = length(t_run)
    n_E = length(E_centers)

    if size(Ie_raw, 1) != n_μ * n_z
        error("Inconsistent dimensions between Ie_raw and (μ_lims, h_atm).")
    end

    # Ie_raw is stacked as [z + (μ-1)*n_z, t, E]. Reshape to [nz, n_μ, nt, nE].
    return reshape(Ie_raw, n_z, n_μ, n_t, n_E)
end


"""
    restructure_streams_of_Ie(Ie, θ_lims, new_θ_lims)

Function that merges the streams of `Ie` that are given over `θ_lims` to fit the
`new_θ_lims` of interest. It can be useful when wanting to merge some streams for plotting.

For example, if we have `θ_lims = [180 160 140 120 100 90 80 60 40 20 0]`, and
we want to plot with `new_θ_lims = [(180, 160), (160, 120)]`, the function will keep the
first stream as is and merge the streams (160°-140°) and (140°-120°) together into a new
stream with limits (160°-120°).

*Important*: The limits in `new_θ_lims` need to match some existing limits in `θ_lims`. In
the example above, `new_θ_lims = [(180, 165)]` would not have worked because 165° is not a
limit that exists in `θ_lims`.

Entries in `new_θ_lims` can be `nothing` to leave a panel empty.

# Calling
`Ie_plot = restructure_streams_of_Ie(Ie, θ_lims, new_θ_lims)`

# Arguments
- `Ie`: array of electron flux with pitch-angle limits `θ_lims`. Of shape [n\\_z, n\\_μ, n\\_t, n\\_E].
- `θ_lims`: pitch-angle limits. Usually a vector or range.
- `new_θ_lims`: new pitch-angle limits. Given as an array of tuples with angles in the range
                0-180° (where 180° is field-aligned down, 0° is field-aligned up). Use `nothing`
                for empty panels. For example:
```
julia> new_θ_lims = [(180, 170)  (170, 150)  (150, 120)  (120, 100)  (100, 90);  # DOWN
                     (0, 10)     (10, 30)    (30, 60)    (60, 80)    (80, 90)]   # UP
```

# Returns
- `Ie_plot`: array of electron flux with the new pitch-angle limits `new_θ_lims`. Of shape
             [n\\_z, n\\_μ\\_new, n\\_t, n\\_E], where n\\_μ\\_new is the number of streams
             in `new_θ_lims`. The second dimension of `Ie_plot` is sorted such that the
             indices go along the first row of `new_θ_lims`, and then the second row.
             In our example with `new_θ_lims` from above, that would be ``[1 2 3 4 5; 6 7 8 9 10]``.

"""
function restructure_streams_of_Ie(Ie, θ_lims, new_θ_lims)
    # Input Ie has shape [n_z, n_μ, n_t, n_E]
    n_μ_new = length(new_θ_lims)
    n_z = size(Ie, 1)
    n_t = size(Ie, 3)
    n_E = size(Ie, 4)
    Ie_plot = Array{eltype(Ie)}(undef, n_z, n_μ_new, n_t, n_E)
    return restructure_streams_of_Ie!(Ie_plot, Ie, θ_lims, new_θ_lims)
end


"""
    restructure_streams_of_Ie!(Ie_plot, Ie, θ_lims, new_θ_lims)

Merge the streams of `Ie` into the preallocated `Ie_plot`, overwriting its contents.
See `restructure_streams_of_Ie` for the pitch-angle bin semantics.
"""
function restructure_streams_of_Ie!(Ie_plot, Ie, θ_lims, new_θ_lims)
    # Input Ie has shape [n_z, n_μ, n_t, n_E]
    n_μ_new = length(new_θ_lims)
    n_z = size(Ie, 1)
    n_t = size(Ie, 3)
    n_E = size(Ie, 4)

    if size(Ie_plot) != (n_z, n_μ_new, n_t, n_E)
        error("`Ie_plot` has incompatible size. Expected $((n_z, n_μ_new, n_t, n_E)), got $(size(Ie_plot)).")
    end

    fill!(Ie_plot, zero(eltype(Ie_plot)))

    # Flatten new_θ_lims from a 2D array to a 1D vector (row by row)
    new_θ_lims_flat = vec(permutedims(new_θ_lims))
    θ_centers = acosd.(mu_avg(θ_lims))

    # Check if all the limits in new_θ_lims match some existing limits in θ_lims.
    θ_lims_sorted = sort(collect(θ_lims))
    for i in eachindex(new_θ_lims_flat)
        θ_bin = new_θ_lims_flat[i]
        isnothing(θ_bin) && continue  # skip empty panels
        if θ_bin[1] ∉ θ_lims
            error("The limit $(θ_bin[1])° in `angles_to_plot` does not match any limit in the simulation's θ_lims.\n" *
                  "Available limits: $(θ_lims_sorted)")
        elseif θ_bin[2] ∉ θ_lims
            error("The limit $(θ_bin[2])° in `angles_to_plot` does not match any limit in the simulation's θ_lims.\n" *
                  "Available limits: $(θ_lims_sorted)")
        end
    end

    # Precompute which original streams contribute to each requested output stream.
    idx_streams = Vector{Vector{Int}}(undef, n_μ_new)
    for i in eachindex(new_θ_lims_flat)
        θ_bin = new_θ_lims_flat[i]
        if isnothing(θ_bin)
            idx_streams[i] = Int[]
            continue
        end

        θ_min = minimum(θ_bin)
        θ_max = maximum(θ_bin)
        idx_streams[i] = findall(θ -> θ_min <= θ <= θ_max, θ_centers)
    end

    # Restructure to [n_z, n_μ_new, n_t, n_E]
    @views for i in eachindex(idx_streams)
        isempty(idx_streams[i]) && continue
        for j in idx_streams[i]
            Ie_plot[:, i, :, :] .+= Ie[:, j, :, :]
        end
    end

    return Ie_plot
end


"""
    interpolate_profile(data_values, data_altitude_km, target_altitude_m;
                       log_interpolation=true)

Interpolate a single profile from one altitude grid to another.

This is a helper function for interpolating individual data profiles.
Interpolation can be performed in linear or logarithmic space.

# Arguments
- `data_values::Vector`: The data to interpolate (e.g., density or temperature)
- `data_altitude_km::Vector`: Altitude grid of the input data (km)
- `target_altitude_m::Vector`: Target altitude grid for interpolation (m)

# Keyword Arguments
- `log_interpolation::Bool=true`: If `true`, interpolation is done in log space (exponential
    extrapolation). Recommended for densities. Use `false` for temperatures.

# Returns
- `interpolated::Vector`: Interpolated data on the target altitude grid
"""
function interpolate_profile(data_values, data_altitude_km, target_altitude_m;
                             log_interpolation = true)
    # Convert target altitude from meters to kilometers
    target_altitude_km = target_altitude_m / 1e3

    # Prepare data for interpolation
    if log_interpolation
        # Use log interpolation for exponential-like profiles (densities)
        interpolator = PCHIPInterpolation(log.(data_values), data_altitude_km;
                                          extrapolation = ExtrapolationType.Linear)
        interpolated = exp.(interpolator(target_altitude_km))
    else
        # Use linear interpolation (e.g., for temperatures)
        interpolator = PCHIPInterpolation(data_values, data_altitude_km;
                                          extrapolation = ExtrapolationType.Linear)
        interpolated = interpolator(target_altitude_km)
    end

    # Check for negative values and throw error if found
    if any(interpolated .< 0)
        error("Interpolation resulted in negative values. Check input data or interpolation method.")
    end

    return interpolated
end
