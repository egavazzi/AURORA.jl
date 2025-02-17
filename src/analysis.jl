using MAT
using LoopVectorization

"""
    make_density_file(directory_to_process)

This function reads into a folder `directory_to_process` containing results from an AURORA.jl
simulation. It loads the particle flux `Ie` (#e⁻/m²/s), calculates the superthermal e-
density `n_e` (#e⁻/m³) from it, and saves `n_e` into a new file "superthermal\\_e\\_density.mat".

The particle flux `Ie` is defined along a magnetic field line and over an (Energy, pitch\\_angle)-grid.
The number density `n_e` calculated is given along a magnetic field line and over an energy
grid. That way, we have the density of electrons with a certain energy at a specific
altitude and time.

# Calling
`make_density_file(directory_to_process)`

# Inputs
- `directory_to_process`: relative path to the simulation folder to process. Example: "Visions2/Alfven_536s"
"""
function make_density_file(directory_to_process)
    println("Calculating the densities from integrating Ie.")
    ## Find the files to process
    full_path_to_directory = pkgdir(AURORA, "data", directory_to_process)
    files = readdir(full_path_to_directory, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]

    ## Extract simulation grid
    f = matopen(files_to_process[1])
        E = read(f, "E")
        μ_lims = vec(read(f, "mu_lims"))
        t_run = read(f, "t_run")
        h_atm = read(f, "h_atm")
    close(f)
    ## Calculate E_middle_bin from E-grid
    E_grid = E
    ΔE = diff(E_grid); ΔE = [ΔE; ΔE[end]]
    E_middle_bin = E_grid + ΔE / 2
    ## Calculate the velocities
    v = v_of_E.(E_middle_bin)

    # Initiate variables
    n_e_local = zeros(length(h_atm), length(t_run), length(E_middle_bin));
    n_e = zeros(length(h_atm), length(t_run), length(E_middle_bin));
    t = t_run
    ## Calculate the densities
    for (i_file, file) in enumerate(files_to_process)
        if i_file > 1
            fill!(n_e_local, 0.0)
        end

        f = matopen(file)
            Ie = read(f, "Ie_ztE") # [n_μ * nz, nt, nE]
            t_run = read(f, "t_run")
        close(f)

        # calculate_density_from_Ie!(h_atm, t_run, μ_lims, E_middle_bin, v, Ie, n_e_local) # when Ie is in #e-/m²/s
        calculate_density_from_Ie!(h_atm, t_run, μ_lims, E_middle_bin, v, Ie .* 1e4, n_e_local) # when Ie is in #e-/cm²/s

        if i_file == 1
            n_e .= n_e_local
            t .= t_run
        elseif i_file > 1
            n_e = hcat(n_e, n_e_local[:, 2:end, :])
            t = vcat(t, t_run[2:end])
        end
        println("File $i_file/", length(files_to_process), " done.")
    end

    ## And save the densities into a file
    savefile = joinpath(full_path_to_directory, "superthermal_e_density.mat")
    f = matopen(savefile, "w")
        write(f, "n_e", n_e)
        write(f, "E", E)
        write(f, "t", t)
        write(f, "h_atm", h_atm)
    close(f)
end

"""
    calculate_density_from_Ie!(h_atm, t_run, μ_lims, E_middle_bin, v, Ie, n_e)

This function converts a particle flux `Ie` (#e⁻/m²/s) into a number density `n_e` (#e⁻/m³).

The particle flux `Ie` is defined along a magnetic field line and over an (Energy, pitch\\_angle)-grid.
The number density `n_e` calculated is given along a magnetic field line and over an energy
grid. That way, we have the density of electrons with a certain energy at a specific
altitude and time.

# Calling
`calculate_density_from_Ie!(h_atm, t_run, μ_lims, E_middle_bin, v, Ie, n_e)`

# Inputs
- `h_atm`: altitude (m), vector [nz]
- `t_run`: time (s), vector [nt]
- `μ_lims`: cosine of the pitch angle limits of the e- beams, vector [n_beam + 1]
- `E_middle_bin`: middle energy of the energy bins (eV), vector [nE]
- `v`: velocity corresponding to the `E_middle_bin` (m/s), vector [nE]
- `Ie`: electron flux (#e⁻/m²/s), 3D array [n_beam * nz, nt, nE]
- `n_e`: electron density (#e⁻/m³), **empty** 3D array [nz, nt, nE]
"""
function calculate_density_from_Ie!(h_atm, t_run, μ_lims, E_middle_bin, v, Ie, n_e)
    Threads.@threads for i_t in eachindex(t_run)
        for i_μ in 1:(length(μ_lims) - 1)
            for i_z in eachindex(h_atm)
                for i_E in eachindex(E_middle_bin)
                    n_e[i_z, i_t, i_E] += 1 / v[i_E] *
                                            Ie[(i_μ - 1) * length(h_atm) + i_z, i_t, i_E]
                end
            end
        end
    end
end


## ====================================================================================== ##


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
- `directory_to_process`: absolute path to the directory to process.
- `downsampling_factor`: downsampling factor for the time

# Outputs
The downsampled electron fluxes `Ie` will be saved in a subfolder inside the `directory_to_process`.
"""
function downsampling_fluxes(directory_to_process, downsampling_factor)
    files = readdir(directory_to_process, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]

    for j in files_to_process
        f = matopen(j)
            Ie = read(f, "Ie_ztE") # [n_μ * nz, nt, nE]
            E = read(f, "E")
            t_run = read(f, "t_run")
            μ_lims = read(f, "mu_lims")
            h_atm = read(f, "h_atm")
            I0 = read(f, "I0")
            μ_scatterings = read(f, "mu_scatterings")
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
            write(file, "E", E)
            write(file, "t_run", t_run)
            write(file, "mu_lims", μ_lims)
            write(file, "h_atm", h_atm)
            write(file, "I0", I0)
            write(file, "mu_scatterings", μ_scatterings)
        close(file)
    end

    return nothing
end


## ====================================================================================== ##


"""
    make_volume_excitation_file(directory_to_process)

Reads into a folder `directory_to_process` containing results from an AURORA.jl simulation,
loads the particle flux `Ie` (#e⁻/m²/s), and calculates the volume-excitation-rates. For
prompt emissions, volume-excitation-rates correspond also to volume-emission-rates.

# Calling
`make_volume_excitation_file(directory_to_process)`

# Inputs
- `directory_to_process`: Name of the simulation folder to process.
    Must be situated under "data/". Example: "Visions2/Alfven_536s"
"""
function make_volume_excitation_file(directory_to_process)
    ## Find the files to process
    full_path_to_directory = pkgdir(AURORA, "data", directory_to_process)
    files = readdir(full_path_to_directory, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]

    ## Load simulation grid
    f = matopen(files_to_process[1])
        h_atm = read(f, "h_atm")
        E = read(f, "E")
    close(f)
    dE = diff(E); dE = [dE; dE[end]]

    ## Load simulation neutral densities
    data = matread(joinpath(full_path_to_directory, "neutral_atm.mat"))
    nO = data["nO"]
    nO2 = data["nO2"]
    nN2 = data["nN2"]

    ## Load the emission cross-sections
    σ_4278 = excitation_4278(E .+ dE/2)
    σ_6730 = excitation_6730_N2(E .+ dE/2)
    σ_7774_O = excitation_7774_O(E .+ dE/2)
    σ_7774_O2 = excitation_7774_O2(E .+ dE/2)
    σ_8446_O = excitation_8446_O(E .+ dE/2)
    σ_8446_O2 = excitation_8446_O2(E .+ dE/2)
    σ_O1D = excitation_O1D(E .+ dE/2)
    σ_O1S = excitation_O1S(E .+ dE/2)
    ## Load/calculate the ionization cross-sections
    σ_N2, σ_O2, σ_O = load_cross_sections(E, dE) # load the cross-sections
    N2_levels, O2_levels, O_levels = load_excitation_threshold() # load the energy levels
    σ_Oi = σ_O' * O_levels[:, 2] # basically do sum(cross_section_for_each_reaction * number_of_ionizations_per_reaction)
    σ_O2i = σ_O2' * O2_levels[:, 2]
    σ_N2i = σ_N2' * N2_levels[:, 2]

    ## Initialize arrays to store the results for each time-slice
    Q4278 = Vector{Matrix{Float64}}()
    Q6730 = Vector{Matrix{Float64}}()
    Q7774_O = Vector{Matrix{Float64}}()
    Q7774_O2 = Vector{Matrix{Float64}}()
    Q7774 = Vector{Matrix{Float64}}()
    Q8446_O = Vector{Matrix{Float64}}()
    Q8446_O2 = Vector{Matrix{Float64}}()
    Q8446 = Vector{Matrix{Float64}}()
    QO1D = Vector{Matrix{Float64}}()
    QO1S = Vector{Matrix{Float64}}()
    QOi = Vector{Matrix{Float64}}()
    QO2i = Vector{Matrix{Float64}}()
    QN2i = Vector{Matrix{Float64}}()
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
            Ie_ztE = read(f, "Ie_ztE")[:, 2:end, :]
            t_local = read(f, "t_run")[2:end]
            close(f)
        end

        ## Sum Ie over the beams
        #=
        I had the idea of initializing Ie_ztE_omni outside of the file-reading loop and then only
        fill it with zeros at the beginning of each new loop, to minimize allocations and GC.
        But then I realized that the current version of the code with the initialization inside
        the loop allows for the Ie_ztE slices to have different n_t sizes.
        Currently (v0.4.2), the Ie_ztE saved in each IeFlickering-XX.mat file all have the
        same size so it is not of big use. But this is not a bottleneck anyway and it allows
        for using time slices with different n_t sizes in the future. // EG 20241221
        =#
        n_z = length(h_atm)
        n_μ = size(Ie_ztE, 1) ÷ n_z # (÷ returns an Int)
        n_t = size(Ie_ztE, 2)
        n_E = size(Ie_ztE, 3)
        Ie_ztE_omni = zeros(n_z, n_t, n_E)
        @views for i_μ in 1:n_μ
            Ie_ztE_omni .+= Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :] # bottleneck, ~75% of time spent here
        end

        # This was an attempt to make things faster with @turbo, but I didn't observe big
        # differences when benchmarking on my machine.
        # @turbo for i_E in 1:n_E
        #     for i_t in 1:n_t
        #         for i_μ in 1:n_μ
        #             for i_z in 1:n_z
        #                 Ie_ztE_omni[i_z, i_t, i_E] += Ie_ztE[(i_μ - 1) * n_z + i_z, i_t, i_E]
        #             end
        #         end
        #     end
        # end

        ## Calculate Q (volume-excitation-rate) for various optical emissions
        Q4278_local = calculate_volume_excitation(h_atm, t_local, Ie_ztE_omni, σ_4278, nN2)
        Q6730_local = calculate_volume_excitation(h_atm, t_local, Ie_ztE_omni, σ_6730, nN2)
        Q7774_O_local = calculate_volume_excitation(h_atm, t_local, Ie_ztE_omni, σ_7774_O, nO)
        Q7774_O2_local = calculate_volume_excitation(h_atm, t_local, Ie_ztE_omni, σ_7774_O2, nO2)
        Q7774_local = Q7774_O_local + Q7774_O2_local
        Q8446_O_local = calculate_volume_excitation(h_atm, t_local, Ie_ztE_omni, σ_8446_O, nO)
        Q8446_O2_local = calculate_volume_excitation(h_atm, t_local, Ie_ztE_omni, σ_8446_O2, nO2)
        Q8446_local = Q8446_O_local + Q8446_O2_local
        QO1D_local = calculate_volume_excitation(h_atm, t_local, Ie_ztE_omni, σ_O1D, nO) # quenching is not taken into account
        QO1S_local = calculate_volume_excitation(h_atm, t_local, Ie_ztE_omni, σ_O1S, nO) # quenching is not taken into account
        # Calculate Q (volume-excitation-rate) for ionizations
        QOi_local = calculate_volume_excitation(h_atm, t_local, Ie_ztE_omni, σ_Oi, nO)
        QO2i_local = calculate_volume_excitation(h_atm, t_local, Ie_ztE_omni, σ_O2i, nO2)
        QN2i_local = calculate_volume_excitation(h_atm, t_local, Ie_ztE_omni, σ_N2i, nN2)

        ## Push the newly calculated Q_local for the current time-slice into a vector
        # We get something like Q4278 = [[n_z, n_t1], [n_z, n_t2], ...]
        push!(Q4278, Q4278_local)
        push!(Q6730, Q6730_local)
        push!(Q7774_O, Q7774_O_local)
        push!(Q7774_O2, Q7774_O2_local)
        push!(Q7774, Q7774_local)
        push!(Q8446_O, Q8446_O_local)
        push!(Q8446_O2, Q8446_O2_local)
        push!(Q8446, Q8446_local)
        push!(QO1D, QO1D_local)
        push!(QO1S, QO1S_local)
        push!(QOi, QOi_local)
        push!(QO2i, QO2i_local)
        push!(QN2i, QN2i_local)
        push!(t, t_local)

        next!(p)
    end

    ## Concatenate along time
    # We get Q4278 = [n_z, n_t]
    Q4278 = reduce(hcat, Q4278)
    Q6730 = reduce(hcat, Q6730)
    Q7774_O = reduce(hcat, Q7774_O)
    Q7774_O2 = reduce(hcat, Q7774_O2)
    Q7774 = reduce(hcat, Q7774)
    Q8446_O = reduce(hcat, Q8446_O)
    Q8446_O2 = reduce(hcat, Q8446_O2)
    Q8446 = reduce(hcat, Q8446)
    QO1D = reduce(hcat, QO1D)
    QO1S = reduce(hcat, QO1S)
    QOi = reduce(hcat, QOi)
    QO2i = reduce(hcat, QO2i)
    QN2i = reduce(hcat, QN2i)
    t = reduce(vcat, t)

    ## Save results
    savefile = joinpath(full_path_to_directory, "Qzt_all_L.mat")
    f = matopen(savefile, "w")
        write(f, "h_atm", h_atm)
        write(f, "t", t)
        write(f, "Q4278", Q4278)
        write(f, "Q6730", Q6730)
        write(f, "Q7774_O", Q7774_O)
        write(f, "Q7774_O2", Q7774_O2)
        write(f, "Q7774", Q7774)
        write(f, "Q8446_O", Q8446_O)
        write(f, "Q8446_O2", Q8446_O2)
        write(f, "Q8446", Q8446)
        write(f, "QO1D", QO1D)
        write(f, "QO1S", QO1S)
        write(f, "QOi", QOi)
        write(f, "QO2i", QO2i)
        write(f, "QN2i", QN2i)
    close(f)

    println("Volume excitation rates saved in $savefile")

    return nothing
end


"""
    calculate_volume_excitation(h_atm, t, Ie_ztE_omni, σ, n)

Calculate the volume-excitation-rate for an excitation of interest, produced by the electron
flux `Ie_ztE_omni` that is summed over the beams (omnidirectional).

The excitation of interest is chosen through the cross-section `σ` given to the function.
Note that the neutral density `n` should match the excitation of interest (e.g. use nN2 when
calculating the volume-excitation-rate of the 4278Å optical emission).

# Calling
`Q = calculate_volume_excitation(h_atm, t, Ie_ztE, σ, n)``

# Inputs
- `h_atm`: altitude (m). Vector [n\\_z]
- `t`: time (s). Vector [n\\_t]
- `Ie_ztE_omni`: omnidirectional electron flux (#e⁻/m²/s). 3D array [n\\_z, n\\_t, n\\_E]
- `σ`: excitation cross-section (m⁻²). Vector [n\\_E]
- `n`: density of exciteable atmospheric specie (m⁻³). Vector [n\\_z]
"""
function calculate_volume_excitation(h_atm, t, Ie_ztE_omni, σ, n)
    # Initialize
    n_z = length(h_atm)
    n_t = length(t)
    Q = zeros(n_z, n_t)
    # Calculate Q for each time step
    @views for i_t in eachindex(t)
        Q[:, i_t] .= (Ie_ztE_omni[:, i_t, :] * σ) .* n
    end

    return Q
end


## ====================================================================================== ##


"""
    make_Ie_top_file(directory_to_process)

Reads into a folder `directory_to_process` containing results from an AURORA.jl simulation
and extracts the particle flux `Ie` (#e⁻/m²/s) at the top of the ionosphere (i.e. at the
max altitude used in the simulation).

# Calling
`make_Ie_top_file(directory_to_process)`

# Inputs
- `directory_to_process`: Name of the simulation folder to process.
    Must be situated under "data/". Example: "Visions2/Alfven_536s"
"""
function make_Ie_top_file(directory_to_process)
    ## Find the files to process
    full_path_to_directory = pkgdir(AURORA, "data", directory_to_process)
    files = readdir(full_path_to_directory, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]

    ## Load simulation grid
    f = matopen(files_to_process[1])
        h_atm = read(f, "h_atm")
        E = read(f, "E")
        μ_scatterings = read(f, "mu_scatterings")
    close(f)
    n_z = length(h_atm)
    dE = diff(E); dE = [dE; dE[end]]
    BeamWeight = μ_scatterings["BeamWeight"]

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
            Ie_ztE = read(f, "Ie_ztE")[:, 2:end, :]
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
    Ie_top = Ie_top ./ reshape(dE, (1, 1, :)) ./ BeamWeight # in #e-/m²/s/eV/ster

    ## Save results
    savefile = joinpath(full_path_to_directory, "Ie_top.mat")
    f = matopen(savefile, "w")
        write(f, "E", E)
        write(f, "dE", dE)
        write(f, "BeamW", BeamWeight)
        write(f, "t", t)
        write(f, "Ie_top_raw", Ie_top_raw)
        write(f, "Ie_top", Ie_top)
    close(f)

    println("Top flux saved in $savefile")

    return nothing
end


## ====================================================================================== ##


"""
    make_current_file(directory_to_process)

Reads into a folder `directory_to_process` containing results from an AURORA.jl simulation,
loads the particle flux `Ie` (#e⁻/m²/s) and calculates the field-aligned current-density
and field-aligned energy-flux for each height and through time.

The following variables are saved to a file named *J.mat*:
- `J_up`: Field-aligned current-density in the upward direction. 2D array [n\\_z, n\\_t]
- `J_down`: Field-aligned current-density in the downward direction. 2D array [n\\_z, n\\_t]
- `E_up`: Field-aligned energy-flux (eV/m²/s) in the upward direction. 2D array [n\\_z, n\\_t]
- `E_down`: Field-aligned energy-flux (eV/m²/s) in the downward direction. 2D array [n\\_z, n\\_t]

# Calling
`make_current_file(directory_to_process)`

# Inputs
- `directory_to_process`: Name of the simulation folder to process.
    Must be situated under "data/". Example: "Visions2/Alfven_536s"
"""
function make_current_file(directory_to_process)
    ## Find the files to process
    full_path_to_directory = pkgdir(AURORA, "data", directory_to_process)
    files = readdir(full_path_to_directory, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]

    ## Load simulation grid
    f = matopen(files_to_process[1])
        h_atm = read(f, "h_atm")
        E = read(f, "E")
        μ_lims = read(f, "mu_lims")
    close(f)
    dE = diff(E); dE = [dE; dE[end]]
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
            Ie_ztE = read(f, "Ie_ztE")[:, 2:end, :]
            t_local = read(f, "t_run")[2:end]
            close(f)
        end

        ## Calculate the field aligned currents and energy flux
        n_z = length(h_atm)
        n_μ = length(μ_center)
        n_t = size(Ie_ztE, 2)
        J_up_local = zeros(n_z, n_t)
        J_down_local = zeros(n_z, n_t)
        IeE_up_local = zeros(n_z, n_t)
        IeE_down_local = zeros(n_z, n_t)
        @views for i_μ in 1:n_μ
            if μ_center[i_μ] > 0
                J_up_local .+= q_e * abs(μ_center[i_μ]) .* sum(Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :], dims=3)
                IeE_up_local .+= abs(μ_center[i_μ]) .* sum(Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :] .* reshape(E, (1, 1, :)), dims=3)
            else
                J_down_local .+= q_e * abs(μ_center[i_μ]) .* sum(Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :], dims=3)
                IeE_down_local .+= abs(μ_center[i_μ]) .* sum(Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :] .* reshape(E, (1, 1, :)), dims=3)
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
    savefile = joinpath(full_path_to_directory, "J.mat")
    f = matopen(savefile, "w")
        write(f, "h_atm", h_atm)
        write(f, "t", t)
        write(f, "J_up", J_up)
        write(f, "J_down", J_down)
        write(f, "IeE_up", IeE_up)
        write(f, "IeE_down", IeE_down)
    close(f)

    println("Currents saved in $savefile")

    return nothing
end


## ====================================================================================== ##



"""
    make_column_excitation_file(directory_to_process)

Reads into a folder `directory_to_process` containing results from an AURORA.jl simulation,
loads the volume excitation rates `Q_XXXX` (#excitation/m³/s) contained in the file `Qzt_all_L.mat`
and integrate them in height, taking into account the finite speed of light.

The calculated colum-integrated excitation rates are saved to a file named *I_lambda_of_t.mat*.
The column-integrated excitation rates are named "I_4278, I_6730, ...". They are all vectors
in time (length n\\_t), and have units of (#excitation/m²/s).

Note that the function `make_volume_excitation_file()` needs to be run before this one, as
we need the file `Qzt_all_L.mat` with the volume excitation rates.

# Calling
`make_current_file(directory_to_process)`

# Inputs
- `directory_to_process`: Name of the simulation folder to process.
    Must be situated under "data/". Example: "Visions2/Alfven_536s"
"""
function make_column_excitation_file(directory_to_process)
    ## Find the files to process
    full_path_to_directory = pkgdir(AURORA, "data", directory_to_process)
    Q_file = joinpath(full_path_to_directory, "Qzt_all_L.mat")

    ## Load volume-excitation-rates
    data = matread(Q_file)
    h_atm = data["h_atm"]
    t = data["t"]
    Q4278 = data["Q4278"]
    Q6730 = data["Q6730"]
    Q7774 = data["Q7774"]
    Q7774_O = data["Q7774_O"]
    Q7774_O2 = data["Q7774_O2"]
    Q8446 = data["Q8446"]
    Q8446_O = data["Q8446_O"]
    Q8446_O2 = data["Q8446_O2"]
    QO1D = data["QO1D"]
    QO1S = data["QO1S"]

    ## Integrate in altitude, taking into account finite photon velocity and path length.
    I_4278 = q2colem(t, h_atm, Q4278)
    I_6730 = q2colem(t, h_atm, Q6730)
    I_7774 = q2colem(t, h_atm, Q7774)
    I_7774_O = q2colem(t, h_atm, Q7774_O)
    I_7774_O2 = q2colem(t, h_atm, Q7774_O2)
    I_8446 = q2colem(t, h_atm, Q8446)
    I_8446_O = q2colem(t, h_atm, Q8446_O)
    I_8446_O2 = q2colem(t, h_atm, Q8446_O2)
    I_O1D = q2colem(t, h_atm, QO1D)
    I_O1S = q2colem(t, h_atm, QO1S)

    ## Save results
    savefile = joinpath(full_path_to_directory, "I_lambda_of_t.mat")
    f = matopen(savefile, "w")
        write(f, "t", t)
        write(f, "I_4278", I_4278)
        write(f, "I_6730", I_6730)
        write(f, "I_7774", I_7774)
        write(f, "I_7774_O", I_7774_O)
        write(f, "I_7774_O2", I_7774_O2)
        write(f, "I_8446", I_8446)
        write(f, "I_8446_O", I_8446_O)
        write(f, "I_8446_O2", I_8446_O2)
        write(f, "I_O1D", I_O1D)
        write(f, "I_O1S", I_O1S)
    close(f)

    println("Column excitation rates saved in $savefile")

    return nothing
end

import Interpolations # use import here to avoid conflicts with DataInterpolations.jl
using Integrals
"""
    q2colem(t, h_atm, Q, A = 1, τ = ones(length(h_atm)))

Integrate the volume-excitation-rate (#exc/m³/s) to column-excitation-rate (#exc/m²/s).
Takes into account the time-delay between light emitted at different altitudes. Photons
emitted at at altitude of 200km will arrive at the detector 100e3/3e8 = 0.333 ms later than
electrons emitted at an altitude of 100km. This is a small time-shift, but it is close to
the time-differences corresponding to the phase-shifts between auroral emissions varying at
~10Hz.

The einstein coefficient `A` and effective lifetime `τ` are optional (equal to one by default).

# Calling
`I = q2colem(t, h_atm, Q, A, τ)`

# Inputs
- `h_atm`: altitude (m). Vector [n\\_z]
- `t`: time (s). Vector [n\\_t]
- `Q`: volume-excitation-rate (#exc/m³/s) of the wavelength of interest. 2D array [n\\_z, n\\_t]
- `A`: einstein coefficient (s⁻¹). Scalar (Float or Int)
- `τ`: effective lifetime (s). Vector [n\\_z].

# Output
- `I`: integrated column-excitation-rate (#exc/m²/s) of the wavelength of interest. Vector [n\\_t]
"""
function q2colem(t, h_atm, Q, A = 1, τ = ones(length(h_atm)))

    ## Define constant
    c = 2.99792458e8 # speed of light (m/s)

    ## Apply the effective lifetime and the Einstein coefficient
    Q = Q .* τ .* A

    ## Create the 2D interpolator
    # To explain how this function is working, let's follow an example.
    # Suppose we have the following Q matrix:
    #     julia> Q
    #     5×5 Matrix{Float64}:
    #     0.65745   0.0652896  0.313073   0.4075     0.811552
    #     0.780153  0.530831   0.546205   0.575431   0.707561
    #     0.573467  0.555628   0.567592   0.88814    0.728269
    #     0.865137  0.141949   0.0853166  0.0506817  0.65735
    #     0.707699  0.646255   0.959993   0.932128   0.956778
    # given over the following time and height vectors:
    #     julia> t
    #     5-element Vector{Int64}:
    #         1
    #         2
    #         3
    #         4
    #         5
    #     julia> h_atm
    #     5-element Vector{Int64}:
    #         1
    #         2
    #         3
    #         4
    #         5
    nodes = (h_atm, t)
    itp = Interpolations.interpolate(nodes, Q, Gridded(Linear()))
    itp = Interpolations.extrapolate(itp, 0.0) # allows for extrapolation but set those values to 0.0

    ## Shift data in time
    # We shift the values in time.
    # For each height (index i) and time (index j) position in the Q matrix, we calculate at
    # what time the corresponding photons will reach the bottom of the height column. That new
    # time is given by the formula `t[j] - (h_atm[i] - h_atm[1]) / c`. We use that new time as
    # input to the interpolator `itp`,  which has for effect to "shift" that value in time.
    # Continuing our example from above, and taking c = 1 for simplicity, we get
    #     julia> I = [itp(h_atm[i], (t[j] - (h_atm[i] - h_atm[1]) / 1)) for i in eachindex(h_atm), j in eachindex(t)]
    #     5×5 Matrix{Float64}:
    #     0.65745  0.0652896  0.313073  0.4075    0.811552
    #     0.0      0.780153   0.530831  0.546205  0.575431
    #     0.0      0.0        0.573467  0.555628  0.567592
    #     0.0      0.0        0.0       0.865137  0.141949
    #     0.0      0.0        0.0       0.0       0.707699
    # We can see that the values have been shifted/delayed in time.
    # For the first time t = 1, only photons from the height h_atm = 1 are arriving. For the
    # time t = 2, photons from the height h_atm = 2 that were emitted at time t = 1 are now
    # arriving. Etc. etc.

    # We used a 2D interpolation. As the height are unchanged, we could have also done a 1D
    # interpolation for each height. We use the 2D interpolation for convenience as it allows
    # us to have only one interpolatior object (instead of one per height). This is also
    # closest to the method used in the legacy Matlab code from which this function is
    # inspired.
    I = [itp(h_atm[i], (t[j] - (h_atm[i] - h_atm[1]) / c)) for i in eachindex(h_atm), j in eachindex(t)]

    ## Integrate in height
    # Now that we have the matrix I which contains the information about the number of
    # photons arriving at the bottom of the column for each height and time steps, we can
    # integrate it along the height.
    problem = SampledIntegralProblem(I, h_atm; dim=1)
    method = TrapezoidalRule()
    I_lambda = solve(problem, method)

    return I_lambda.u
end
