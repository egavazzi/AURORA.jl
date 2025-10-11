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
	mₑ = 9.10939e-31;
	qₑ = 1.6021773e-19;

	v = (2 * qₑ * abs(E) / mₑ) .^ (1/2) .* sign(E);
	return v
end

function CFL_criteria(t, h_atm, v, CFL_number=64)
    # The Courant-Freidrichs-Lewy (CFL) number normally hase to be small (<4) to ensure numerical
    # stability. However, as a Crank-Nicolson scheme is always stable, we can take a bigger CFL. We
    # should be careful about numerical accuracy though.
    # For Gaussian inputs (or similar), it seems that the CFL can be set to 64 without major effects
    # on the results, while reducing computational time tremendously
    dt = t[2] - t[1]
    dz = h_atm[2] - h_atm[1]
    # Calculate the maximum dt that still satisfies the CFL criteria.
    dt_max = CFL_number * dz / v
    # Then find the smallest integer that dt can be divided to become smaller than dt_max.
    CFL_factor = ceil(Int, dt / dt_max) # that integer is the CFL_factor
    # Now make t_finer using the CFL_factor
    t_finer = range(t[1], t[end], CFL_factor * (length(t) - 1) + 1)

    return t_finer, CFL_factor
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


## ====================================================================================== ##

import LibGit2
import Pkg
function save_parameters(altitude_lims, θ_lims, E_max, B_angle_to_zenith, t_sampling, t,
    n_loop, CFL_number, INPUT_OPTIONS, savedir)
	savefile = joinpath(savedir, "parameters.txt")
    commit_hash = if isdir(joinpath(pkgdir(AURORA), ".git"))
        LibGit2.head(pkgdir(AURORA))
    else
        "Not available"
    end
    version_AURORA = pkgversion(AURORA)
    open(savefile, "w") do f
        write(f, "altitude_lims = $altitude_lims \n")
        write(f, "θ_lims = $θ_lims \n")
        write(f, "E_max = $E_max \n")
        write(f, "B_angle_to_zenith = $B_angle_to_zenith \n")
        write(f, "\n")
        write(f, "t_sampling = $t_sampling \n")
        write(f, "t = $t \n")
        write(f, "n_loop = $n_loop \n")
        write(f, "\n")
        write(f, "CFL_number = $CFL_number")
        write(f, "\n")
        write(f, "input_options = $INPUT_OPTIONS \n")
        write(f, "\n")
        write(f, "commit_hash = $commit_hash \n")
        write(f, "version_AURORA = $version_AURORA")
    end
end


using MAT: matopen
function save_neutrals(h_atm, n_neutrals, ne, Te, Tn, savedir)
    savefile = joinpath(savedir, "neutral_atm.mat")
    file = matopen(savefile, "w")
        write(file, "h_atm", h_atm)
        write(file, "nN2", n_neutrals.nN2)
        write(file, "nO2", n_neutrals.nO2)
        write(file, "nO", n_neutrals.nO)
        write(file, "ne", ne)
        write(file, "Te", Te)
        write(file, "Tn", Tn)
    close(file)
end


using MAT: matopen
using Printf: @sprintf
function save_results(Ie, E, t, μ_lims, h_atm, I0, μ_scatterings, i, CFL_factor, savedir)
    # Extract the time array for the current loop
	t_run = collect(t .+ t[end] * (i - 1))

    # Reduce t_run and Ie to match the t_sampling
    t_run = t_run[1:CFL_factor:end]
    Ie_save = Ie[:, 1:CFL_factor:end, :]

    savefile = joinpath(savedir, (@sprintf "IeFlickering-%02d.mat" i))
	file = matopen(savefile, "w")
		write(file, "Ie_ztE", Ie_save)
		write(file, "E", E)
		write(file, "t_run", t_run)
		write(file, "mu_lims", μ_lims)
		write(file, "h_atm", h_atm)
		write(file, "I0", I0)
		write(file, "mu_scatterings", μ_scatterings)
	close(file)
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
    # Remove trailing slash (if any)
    savefile_clean = rstrip(savefile, '/')

    # Check if path exists (either as file or directory)
    if !ispath(savefile_clean)
        return savefile_clean
    end

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
    find_Ietop_file(path_to_directory)

Look for Ie\\_incoming file present in the directory given by `path_to_directory`. If several
files are starting with the name "Ie\\_incoming", return an error. If only one file is found,
return a string with the path to that file.

# Calling
`Ietop_file = find_Ietop_file(path_to_directory)`

# Inputs
- `path_to_directory`: path to a directory

# Returns
- `Ietop_file`: path to the Ie\\_incoming file, in the form "path_to_directory/Ie_incoming_*.mat"
"""
function find_Ietop_file(path_to_directory)
    incoming_files = filter(file -> startswith(file, "Ie_incoming_"), readdir(path_to_directory))
    if length(incoming_files) > 1
        error("More than one file contains incoming flux. This is not normal")
    else
        return Ietop_file = joinpath(path_to_directory, incoming_files[1])
    end
end


"""
    make_savedir(root_savedir, name_savedir; behavior = "default")

Return the path to the directory where the results will be saved. If the directory does not
already exist, create it.

If the constructed `savedir` already exists and contains files starting with `"IeFlickering-"`,
a new directory is created to avoid accidental overwriting of results (e.g., `savedir(1)`, `savedir(2)`, etc.).

# Calling
`savedir = make_savedir(root_savedir, name_savedir)` \\
`savedir = make_savedir(root_savedir, name_savedir; behavior = "custom")`

# Arguments
- `root_savedir::String`: The root directory where the data will be saved. If empty or
    contains only spaces, it defaults to `"backup"`.
- `name_savedir::String`: The name of the subdirectory to be created within `root_savedir`.
    If empty or contains only spaces, it defaults to the current date and time in the
    format `"yyyymmdd-HHMM"`.
- `behavior::String` (optional): Determines how the full path is constructed.
    - `"default"`: The path will be built starting under the `data/` folder of the AURORA installation
        (i.e., `AURORA_folder/data/root_savedir/name_savedir/`, where `AURORA_folder` is
        the folder containing the AURORA code). This is the default behavior.
    - `"custom"`: The path will be built as `root_savedir/name_savedir/`, with the argument
        `root_savedir` treated as an absolute or relative path. This allows for saving results
        in any location on the system. Useful if AURORA is installed as a dependency to some
        other project.

# Returns
- `savedir::String`: The full path to the directory where the results will be saved.
"""
function make_savedir(root_savedir, name_savedir; behavior = "default")
    ## Create the folder to save the data to
    # If `root_savedir` is empty or contains only "space" characters, we use "backup/" as a name
    if isempty(root_savedir) || !occursin(r"[^ ]", root_savedir)
        root_savedir = "backup"
    end
    # If `name_savedir` is empty or contains only "space" characters, we use the current date and time as a name
    if isempty(name_savedir) || !occursin(r"[^ ]", name_savedir)
        name_savedir = string(Dates.format(now(), "yyyymmdd-HHMM"))
    end

    # Make a string with full path of savedir from root_savedir and name_savedir
    if behavior == "default"
        savedir = pkgdir(AURORA, "data", root_savedir, name_savedir)
    elseif behavior == "custom"
        savedir = joinpath(root_savedir, name_savedir)
    end

    # Rename `savedir` to `savedir(1)` if it exists and already contain results. If
    # `savedir(1)` exists then it will be renamed to `savedir(2)` and so on
    if isdir(savedir) && (filter(startswith("IeFlickering-"), readdir(savedir)) |> length) > 0
        savedir = rename_if_exists(savedir)
    end

    # And finally create the directory
    mkpath(savedir)
    print("\n", @bold "Results will be saved at $savedir \n")

    return savedir
end

## ====================================================================================== ##


function square(x)
    ifelse(mod2pi(x) < π, 1.0, -1.0)
end


## ====================================================================================== ##

# function from Björn for smooth input onset
function f_smooth_transition(x, a = 0, b = 1)
    if (b - a) == 0
        y = 1
        return y
    end

    x = (x - a) / (b - a)
    psi0P = psi(x)
    psi1N = psi(1 - x)

    y = psi0P / (psi0P + psi1N)

    return y
end

function psi(x)
    PSI = exp(-1 / x)
    if x <= 0
        PSI = 0
    end
    return PSI
end
f_smooth_transition(0, 1, 0)

## ====================================================================================== ##

# Function to restructure the matrix from 3D [n_mu x nz, nt, nE] to 4D [n_mu, nz, nt, nE]
function restructure_Ie_from_3D_to_4D(Ie_raw, μ_lims, h_atm, t_run, E)
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

For example, if we have `θ_lims = [180 160 140 120 100 90 80 60 40 20 0]`, and
we want to plot with `new_θ_lims = [180 160 120 100 90 80 60 40 20 0]` (note that
this should be an array of tuples, but to simplify the comparison we write it as a vector
here), the function will merge the streams (160°-140°) and (140°-120°) together into a new
stream with limits (160°-120°).

*Important*: The limits in `new_θ_lims` need to match some existing limits in `θ_lims`. In
the example above, `new_θ_lims = [180 160 120 100 90 80 65 40 20 0]` would not have worked
because 65° is not a limit that exists in `θ_lims`.

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

    # Check if all the limits in new_θ_lims match some limits in θ_lims
    for i in eachindex(new_θ_lims_temp)
        if new_θ_lims_temp[i][1] ∉ θ_lims
            error("The limit $(new_θ_lims_temp[i][1]) in `new_θ_lims` does not match any limit in `θ_lims`.")
        elseif new_θ_lims_temp[i][2] ∉ θ_lims
            error("The limit $(new_θ_lims_temp[i][2]) in `new_θ_lims` does not match any limit in `θ_lims`.")
        end
    end


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
