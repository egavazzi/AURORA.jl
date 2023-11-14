using MAT

"""
    make_density_file(directory_to_process)

This function reads into a folder `directory_to_process` containing results from an AURORA.jl
simulation. It loads the particle flux `Ie` (#e⁻/cm²/s), calculates the superthermal e-
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

    ## Extracte simulation grid
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

        calculate_density_from_Ie!(h_atm, t_run, μ_lims, E_middle_bin, v, Ie, n_e_local)

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

This function converts a particle flux `Ie` (#e⁻/cm²/s) into a number density `n_e` (#e⁻/m³).

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
- `Ie`: electron flux (#e⁻/cm²/s), 3D array [n_beam * nz, nt, nE]
- `n_e`: electron density (#e⁻/m³), **empty** 3D array [nz, nt, nE]
"""
function calculate_density_from_Ie!(h_atm, t_run, μ_lims, E_middle_bin, v, Ie, n_e)
    Threads.@threads for i_t in eachindex(t_run)
        for i_μ in 1:(length(μ_lims) - 1)
            for i_z in eachindex(h_atm)
                for i_E in eachindex(E_middle_bin)
                    n_e[i_z, i_t, i_E] += 1 / v[i_E] *
                                            Ie[(i_μ - 1) * length(h_atm) + i_z, i_t, i_E] *
                                            1e4 # to convert from cm² to m²
                end
            end
        end
    end
end
