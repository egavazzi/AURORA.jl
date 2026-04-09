using MAT: matopen, matread
using ProgressMeter: Progress, next!


# ======================================================================================== #
#                              FLUX EXTRACTION & DOWNSAMPLING                            #
# ======================================================================================== #

"""
    make_Ie_top_file(directory_to_process)

Reads into a folder `directory_to_process` containing results from an AURORA.jl simulation
and extracts the particle flux `Ie` (#e⁻/m²/s) at the top of the ionosphere (i.e. at the
max altitude used in the simulation).

# Calling
`make_Ie_top_file(directory_to_process)`

# Inputs
- `directory_to_process`: absolute or relative path to the simulation directory to process.
"""
function make_Ie_top_file(directory_to_process)
    ## Find the files to process
    files = readdir(directory_to_process, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]
    # The files are sorted in lexicographical order, so IeFlickering-100.mat will be loaded
    # before "IeFlickering-11.mat. We fix that with the following line which sorts them by
    # the number in the filename.
    sort!(files_to_process, by = x -> parse(Int, match(r"IeFlickering-(\d+)\.mat", basename(x))[1]))

    ## Load simulation grid
    f = matopen(files_to_process[1])
        z = read(f, "h_atm")
        E_centers = read(f, "E_centers")
        ΔE = read(f, "dE")
        scattering_data = read(f, "mu_scatterings")
    close(f)
    n_z = length(z)
    Ω_beam = scattering_data["BeamWeight"]

    ## Initialize arrays to store the results for each time-slice
    Ie_top = Vector{Array{Float64, 3}}()
    t = Vector{}()

    n_files = length(files_to_process)
    p = Progress(n_files; desc=string("Processing data"), dt=1.0, color=:blue)
    ## Loop over the files
    for (i_file, file) in enumerate(files_to_process)
        ## Load simulation results for current file.
        if i_file == 1
            f = matopen(file)
                Ie_ztE = read(f, "Ie_ztE")
                t_local = read(f, "t_run")
            close(f)
        else
            f = matopen(file)
            @views Ie_ztE = read(f, "Ie_ztE")[:, 2:end, :]
            t_local = read(f, "t_run")[2:end]
            close(f)
        end

        ## Extract Ie at the top for each beam
        Ie_top_local = Ie_ztE[n_z:n_z:end, :, :]

        ## Push the Ie_top of the current time-slice into a vector
        # We get Ie_top = [[n_μ, n_t1, n_E], [n_μ, n_t2, n_E], ...]
        push!(Ie_top, Ie_top_local)
        push!(t, t_local)

        next!(p)
    end

    ## Concatenate along time
    # We get Ie_top = [n_μ, n_t, n_E]
    Ie_top = reduce(hcat, Ie_top)
    t = reduce(vcat, t)

    ## Play with the units
    Ie_top_raw = copy(Ie_top) # in #e-/m²/s
    Ie_top = Ie_top ./ reshape(ΔE, (1, 1, :)) ./ Ω_beam # in #e-/m²/s/eV/ster

    ## Save results
    savefile = joinpath(directory_to_process, "Ie_top.mat")
    f = matopen(savefile, "w")
        write(f, "E_centers", E_centers)
        write(f, "dE", ΔE)
        write(f, "BeamW", Ω_beam)
        write(f, "t", t)
        write(f, "Ie_top_raw", Ie_top_raw)
        write(f, "Ie_top", Ie_top)
    close(f)

    println("Top flux saved in $savefile")

    return nothing
end


"""
    downsampling_fluxes(directory_to_process, downsampling_factor)

This function extracts `Ie` from the simulation results in `directory_to_process` and
downsample it in time.
For example: if `Ie` is given with a time step of 1ms and we use a `downsampling_factor` of
10, this function will extract the values of `Ie` with a time step of 10ms. It will then
save the results in a new subfolder called`downsampled_10x`, inside the `directory_to_process`.

# Calling
`downsampling_fluxes(directory_to_process, downsampling_factor)`

# Inputs
- `directory_to_process`: absolute or relative path to the simulation directory to process.
- `downsampling_factor`: downsampling factor for the time

# Outputs
The downsampled electron fluxes `Ie` will be saved in a subfolder inside the `directory_to_process`.
"""
function downsampling_fluxes(directory_to_process, downsampling_factor)
    files = readdir(directory_to_process, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]
    # The files are sorted in lexicographical order, so IeFlickering-100.mat will be loaded
    # before "IeFlickering-11.mat. We fix that with the following line which sorts them by
    # the number in the filename.
    sort!(files_to_process, by = x -> parse(Int, match(r"IeFlickering-(\d+)\.mat", basename(x))[1]))

    for j in files_to_process
        f = matopen(j)
            Ie = read(f, "Ie_ztE") # [n_μ * nz, nt, nE]
            E_centers = read(f, "E_centers")
            E_edges = read(f, "E_edges")
            ΔE = read(f, "dE")
            t_run = read(f, "t_run")
            μ_lims = read(f, "mu_lims")
            z = read(f, "h_atm")
            I0 = read(f, "I0")
            scattering_data = read(f, "mu_scatterings")
        close(f)

        # downsample the data
        dt = diff(t_run)[1]
        println("The time-step from simulation is ", dt, "s.")
        new_dt = dt * downsampling_factor
        println("The time-step of the new file will be ", new_dt, "s.")
        Ie = Ie[:, 1:downsampling_factor:end, :]
        t_run = t_run[1:downsampling_factor:end]

        # create new subdir if it doesn't exist
        new_subdir = "downsampled_" * string(downsampling_factor) * "x"
        full_path_to_new_subdir = joinpath(directory_to_process, new_subdir)
        mkpath(full_path_to_new_subdir)

        # create new file
        new_filename = splitdir(j)[2][1:end-4] * "d.mat" # add a 'd' after the number to indicate that it is downsampled
        full_path_to_new_filename = joinpath(full_path_to_new_subdir, new_filename)
        file = matopen(full_path_to_new_filename, "w")
            write(file, "Ie_ztE", Ie)
            write(file, "E_centers", E_centers)
            write(file, "E_edges", E_edges)
            write(file, "dE", ΔE)
            write(file, "t_run", t_run)
            write(file, "mu_lims", μ_lims)
            write(file, "h_atm", z)
            write(file, "I0", I0)
            write(file, "mu_scatterings", scattering_data)
        close(file)
    end

    return nothing
end


# ======================================================================================== #
#                                  FIELD-ALIGNED CURRENTS                                #
# ======================================================================================== #

"""
    make_current_file(directory_to_process)

Reads into a folder `directory_to_process` containing results from an AURORA.jl simulation,
loads the particle flux `Ie` (#e⁻/m²/s) and calculates the field-aligned current-density
and field-aligned energy-flux for each height and through time.

The following variables are saved to a file named *J.mat*:
- `J_up`: Field-aligned current-density in the upward direction. 2D array [n\\_z, n\\_t]
- `J_down`: Field-aligned current-density in the downward direction. 2D array [n\\_z, n\\_t]
- `IeE_up`: Field-aligned energy-flux (eV/m²/s) in the upward direction. 2D array [n\\_z, n\\_t]
- `IeE_down`: Field-aligned energy-flux (eV/m²/s) in the downward direction. 2D array [n\\_z, n\\_t]

# Calling
`make_current_file(directory_to_process)`

# Inputs
- `directory_to_process`: absolute or relative path to the simulation directory to process.
"""
function make_current_file(directory_to_process)
    ## Find the files to process
    files = readdir(directory_to_process, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]
    # The files are sorted in lexicographical order, so IeFlickering-100.mat will be loaded
    # before "IeFlickering-11.mat. We fix that with the following line which sorts them by
    # the number in the filename.
    sort!(files_to_process, by = x -> parse(Int, match(r"IeFlickering-(\d+)\.mat", basename(x))[1]))

    ## Load simulation grid
    f = matopen(files_to_process[1])
        z = read(f, "h_atm")
        E_centers = read(f, "E_centers")
        μ_lims = read(f, "mu_lims")
    close(f)
    μ_center = mu_avg(acosd.(μ_lims))

    ## Initialize arrays to store the results for each time-slice
    J_up = Vector{Array{Float64, 2}}()
    J_down = Vector{Array{Float64, 2}}()
    IeE_up = Vector{Array{Float64, 2}}()
    IeE_down = Vector{Array{Float64, 2}}()
    t = Vector{}()

    ## Define constant
    q_e = 1.602176620898e-19 # elementary charge (C)

    n_files = length(files_to_process)
    p = Progress(n_files; desc=string("Processing data"), dt=1.0, color=:blue)
    ## Loop over the files
    for (i_file, file) in enumerate(files_to_process)
        ## Load simulation results for current file.
        if i_file == 1
            f = matopen(file)
                Ie_ztE = read(f, "Ie_ztE")
                t_local = read(f, "t_run")
            close(f)
        else
            f = matopen(file)
            @views Ie_ztE = read(f, "Ie_ztE")[:, 2:end, :]
            t_local = read(f, "t_run")[2:end]
            close(f)
        end

        ## Calculate the field aligned currents and energy flux
        n_z = length(z)
        n_μ = length(μ_center)
        n_t = size(Ie_ztE, 2)
        J_up_local = zeros(n_z, n_t)
        J_down_local = zeros(n_z, n_t)
        IeE_up_local = zeros(n_z, n_t)
        IeE_down_local = zeros(n_z, n_t)
        @views for i_μ in 1:n_μ
            if μ_center[i_μ] > 0
                J_up_local .+= q_e * abs(μ_center[i_μ]) .* sum(Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :], dims=3)
                IeE_up_local .+= abs(μ_center[i_μ]) .* sum(Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :] .* reshape(E_centers, (1, 1, :)), dims=3)
            else
                J_down_local .+= q_e * abs(μ_center[i_μ]) .* sum(Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :], dims=3)
                IeE_down_local .+= abs(μ_center[i_μ]) .* sum(Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :] .* reshape(E_centers, (1, 1, :)), dims=3)
            end
        end

        ## Push the J of the current time-slice into a vector
        # We get J_up = [[n_z, n_t1], [n_z, n_t2], ...]
        push!(J_up, J_up_local)
        push!(J_down, J_down_local)
        push!(IeE_up, IeE_up_local)
        push!(IeE_down, IeE_down_local)
        push!(t, t_local)

        next!(p)
    end

    ## Concatenate along time
    # We get J_up = [n_z, n_t]
    J_up = reduce(hcat, J_up)
    J_down = reduce(hcat, J_down)
    IeE_up = reduce(hcat, IeE_up)
    IeE_down = reduce(hcat, IeE_down)
    t = reduce(vcat, t)

    ## Save results
    savefile = joinpath(directory_to_process, "J.mat")
    f = matopen(savefile, "w")
        write(f, "h_atm", z)
        write(f, "t", t)
        write(f, "J_up", J_up)
        write(f, "J_down", J_down)
        write(f, "IeE_up", IeE_up)
        write(f, "IeE_down", IeE_down)
    close(f)

    println("Currents saved in $savefile")

    return nothing
end


# ======================================================================================== #
#                         AuroraSimulation convenience wrappers                          #
# ======================================================================================== #

"""
    make_Ie_top_file(sim::AuroraSimulation)

Convenience wrapper that calls [`make_Ie_top_file`](@ref) on `sim.savedir`.
"""
make_Ie_top_file(sim::AuroraSimulation) = make_Ie_top_file(sim.savedir)

"""
    make_current_file(sim::AuroraSimulation)

Convenience wrapper that calls [`make_current_file`](@ref) on `sim.savedir`.
"""
make_current_file(sim::AuroraSimulation) = make_current_file(sim.savedir)

"""
    downsampling_fluxes(sim::AuroraSimulation, downsampling_factor)

Convenience wrapper that calls [`downsampling_fluxes`](@ref) on `sim.savedir`.
"""
downsampling_fluxes(sim::AuroraSimulation, downsampling_factor) = downsampling_fluxes(sim.savedir, downsampling_factor)
