using Dates: Dates, now
using HCubature: hcubature
using Interpolations: linear_interpolation
using MAT: matopen
using ProgressMeter: Progress, next!
using QuadGK: quadgk


# ======================================================================================== #
#                                   UTILITY FUNCTIONS                                     #
# ======================================================================================== #


"""
    find_cascading_file(E_grid, species_dir)

Search for a pre-computed cascading spectra file with matching energy grid.

# Arguments
- `E_grid`: Energy grid to match
- `species_dir`: Directory containing cascading data files for the species

# Returns
- `(file_found, filepath)`: Tuple of boolean and filepath string
"""
function find_cascading_file(E_grid, species_dir)
    cascading_files = readdir(species_dir)

    for filename in cascading_files
        if !isdir(filename)
            try
                filepath = joinpath(species_dir, filename)
                file = matopen(filepath)
                E_grid_saved = read(file, "E4Q")
                close(file)

                # Check if saved grid matches requested grid
                if length(E_grid) <= length(E_grid_saved) && E_grid_saved[1:length(E_grid)] == E_grid
                    return (true, filepath)
                end
            catch
                continue
            end
        end
    end

    return (false, "")
end


"""
    load_cascading_matrices(filepath)

Load pre-computed cascading matrices from a file.

# Arguments
- `filepath`: Path to the .mat file containing cascading data

# Returns
- `(Q_transfer_matrix, E_grid_for_Q, ionization_thresholds)`: Tuple of loaded data
"""
function load_cascading_matrices(filepath)
    println("Loading cascading-matrices from file: ", basename(filepath))
    file = matopen(filepath)
    Q_transfer_matrix = read(file, "Q")
    ionization_thresholds = read(file, "E_ionizations")
    E_grid_for_Q = read(file, "E4Q")
    close(file)

    return (Q_transfer_matrix, E_grid_for_Q, ionization_thresholds)
end


"""
    save_cascading_matrices(Q_transfer_matrix, E_grid_for_Q, ionization_thresholds, species_dir, species_name)

Save calculated cascading matrices to a file.

# Arguments
- `Q_transfer_matrix`: Transfer matrix to save
- `E_grid_for_Q`: Energy grid used for calculations
- `ionization_thresholds`: Ionization threshold energies
- `species_dir`: Directory to save the file
- `species_name`: Name of the species (for filename)
"""
function save_cascading_matrices(Q_transfer_matrix, E_grid_for_Q, ionization_thresholds, species_dir, species_name)
    filename = joinpath(species_dir,
                       string("CascadingSpec", species_name, "ionization_",
                             Dates.format(now(), "yyyymmdd-HHMMSS"),
                             ".mat"))
    file = matopen(filename, "w")
    write(file, "Q", Q_transfer_matrix)
    write(file, "E4Q", E_grid_for_Q)
    write(file, "E_ionizations", ionization_thresholds)
    close(file)
end










# ======================================================================================== #
#                                  CASCADING FUNCTIONS                                    #
# ======================================================================================== #


let Q_transfer_matrix = [], E_grid_for_Q = [], ionization_thresholds = []
    """
        cascading_N2(E_grid, E_primary_energy, E_ionization_threshold, spectrum_type)

    Compute cascading electron spectra for N₂ ionization.

    # Arguments
    - `E`: Energy grid for electrons (eV), Vector [nE] (last element of the grid is omitted)
    - `dE`: Energy grid step size (eV), Vector [nE]
    - `E_primary_energy`: Primary electron energy (eV), scalar
    - `E_ionization_threshold`: Ionization threshold energy (eV), scalar
    - `spectrum_type`: "s" for secondary electrons, "c" for cascading (degraded primary) electrons

    # Returns
    - Electron energy spectrum, unnormalized, Vector [nE]. Depending on `spectrum_type`, it is
    either the spectra of the secondary e⁻, or of the degraded primary e⁻

    # Reference
    Equation from Itikawa 1986 J. Phys. Chem. Ref. Data
    """
    global function cascading_N2(E_grid, dE, E_primary_energy, E_ionization_threshold, spectrum_type)

        # Width of the Lorentzian like secondary spectra, from Itikawa 1986 (eV)
        lorentzian_width = 11.4

        # Check if we need to recalculate or reload the transfer matrix
        # 1. if Q is empty, definitely yes
        # 2. if the E_grid_for_Q is smaller than E_grid, then yes (but the opposite is fine)
        # 3. if the E_grid_for_Q and E_grid are different, then yes
        needs_recalculation = isempty(Q_transfer_matrix) ||
                              length(E_grid) > length(E_grid_for_Q) ||
                              E_grid_for_Q[1:length(E_grid)] != E_grid

        if needs_recalculation
            species_dir = pkgdir(AURORA, "internal_data", "e_cascading", "N2")

            # Try to find a pre-computed cascading spectra file with matching energy grid
            file_found, filepath = find_cascading_file(E_grid, species_dir)

            if file_found
                Q_transfer_matrix, E_grid_for_Q, ionization_thresholds = load_cascading_matrices(filepath)
            else
                println("Could not find a file with matching energy grid.")
                println("Starting to calculate the requested cascading-matrices.")

                # Calculate transfer matrices
                Q_transfer_matrix, E_grid_for_Q, ionization_thresholds = calculate_cascading_N2(E_grid, dE, lorentzian_width)

                # Save the results for future use
                save_cascading_matrices(Q_transfer_matrix, E_grid_for_Q, ionization_thresholds,
                                       species_dir, "N2")
            end
        end

        # Return the requested spectrum type
        if spectrum_type == "s"
            # Secondary electron spectrum: simple Lorentzian formula
            E_excess = E_primary_energy - E_ionization_threshold
            secondary_spectrum = 1 ./ (lorentzian_width^2 .+ E_grid.^2) .* (E_grid .< E_excess / 2)
            return secondary_spectrum

        elseif spectrum_type == "c"
            # Cascading (degraded primary) spectrum: extract directly from transfer matrix
            i_threshold = findmin(abs.(E_ionization_threshold .- ionization_thresholds))[2]
            i_primary = findmin(abs.(E_grid_for_Q .- E_primary_energy))[2]
            cascading_spectrum = Q_transfer_matrix[i_primary, 1:length(E_grid), i_threshold]
            return cascading_spectrum
        end
    end
end

let Q_transfer_matrix = [], E_grid_for_Q = [], ionization_thresholds = []
    """
        cascading_O2(E_grid, dE, E_primary_energy, E_ionization_threshold, spectrum_type)

    Compute cascading electron spectra for O₂ ionization.

    # Arguments
    - `E_grid`: Energy grid for electrons (eV), Vector [nE] (last element of the grid is omitted)
    - `dE`: Energy grid step size (eV), Vector [nE]
    - `E_primary_energy`: Primary electron energy (eV), scalar
    - `E_ionization_threshold`: Ionization threshold energy (eV), scalar
    - `spectrum_type`: "s" for secondary electrons, "c" for cascading (degraded primary) electrons

    # Returns
    - Electron energy spectrum, unnormalized, Vector [nE]. Depending on `spectrum_type`, it is
    either the spectra of the secondary e⁻, or of the degraded primary e⁻
    """
    global function cascading_O2(E_grid, dE, E_primary_energy, E_ionization_threshold, spectrum_type)

        # Width of the Lorentzian like secondary spectra (eV)
        lorentzian_width = 15.2

        # Check if we need to recalculate or reload the transfer matrix
        needs_recalculation = isempty(Q_transfer_matrix) ||
                              length(E_grid) > length(E_grid_for_Q) ||
                              E_grid_for_Q[1:length(E_grid)] != E_grid

        if needs_recalculation
            species_dir = pkgdir(AURORA, "internal_data", "e_cascading", "O2")

            # Try to find a pre-computed cascading spectra file with matching energy grid
            file_found, filepath = find_cascading_file(E_grid, species_dir)

            if file_found
                Q_transfer_matrix, E_grid_for_Q, ionization_thresholds = load_cascading_matrices(filepath)
            else
                println("Could not find a file with matching energy grid.")
                println("Starting to calculate the requested cascading-matrices.")

                # Calculate transfer matrices
                Q_transfer_matrix, E_grid_for_Q, ionization_thresholds = calculate_cascading_O2(E_grid, dE, lorentzian_width)

                # Save the results for future use
                save_cascading_matrices(Q_transfer_matrix, E_grid_for_Q, ionization_thresholds,
                                       species_dir, "O2")
            end
        end

        if spectrum_type == "s"
            E_excess = E_primary_energy - E_ionization_threshold
            secondary_spectrum = 1 ./ (lorentzian_width^2 .+ E_grid.^2) .* (E_grid .< E_excess / 2)
            return secondary_spectrum

        elseif spectrum_type == "c"
            i_threshold = findmin(abs.(E_ionization_threshold .- ionization_thresholds))[2]
            i_primary = findmin(abs.(E_grid_for_Q .- E_primary_energy))[2]
            cascading_spectrum = Q_transfer_matrix[i_primary, 1:length(E_grid), i_threshold]
            return cascading_spectrum
        end
    end
end


let Q_transfer_matrix = [], E_grid_for_Q = [], ionization_thresholds = []
    """
        cascading_O(E_grid, dE, E_primary_energy, E_ionization_threshold, spectrum_type)

    Compute cascading electron spectra for atomic O ionization.

    # Arguments
    - `E_grid`: Energy grid for electrons (eV), Vector [nE] (last element of the grid is omitted)
    - `dE`: Energy grid step size (eV), Vector [nE]
    - `E_primary_energy`: Primary electron energy (eV), scalar
    - `E_ionization_threshold`: Ionization threshold energy (eV), scalar
    - `spectrum_type`: "s" for secondary electrons, "c" for cascading (degraded primary) electrons

    # Returns
    - Electron energy spectrum, unnormalized, Vector [nE]. Depending on `spectrum_type`, it is
    either the spectra of the secondary e⁻, or of the degraded primary e⁻
    """
    global function cascading_O(E_grid, dE, E_primary_energy, E_ionization_threshold, spectrum_type)

        # Interpolate parameters based on primary energy
        A_factor, B_factor = interpolate_O_parameters(E_primary_energy)
        A_factor *= 1.25  # Empirical correction factor

        # Check if we need to recalculate or reload the transfer matrix
        needs_recalculation = isempty(Q_transfer_matrix) ||
                              length(E_grid) > length(E_grid_for_Q) ||
                              E_grid_for_Q[1:length(E_grid)] != E_grid

        if needs_recalculation
            species_dir = pkgdir(AURORA, "internal_data", "e_cascading", "O")

            # Try to find a pre-computed cascading spectra file with matching energy grid
            file_found, filepath = find_cascading_file(E_grid, species_dir)

            if file_found
                Q_transfer_matrix, E_grid_for_Q, ionization_thresholds = load_cascading_matrices(filepath)
            else
                println("Could not find a file with matching energy grid.")
                println("Starting to calculate the requested cascading-matrices.")

                # Calculate transfer matrices
                Q_transfer_matrix, E_grid_for_Q, ionization_thresholds = calculate_cascading_O(E_grid, dE)

                # Save the results for future use
                save_cascading_matrices(Q_transfer_matrix, E_grid_for_Q, ionization_thresholds,
                                       species_dir, "O")
            end
        end

        # Return the requested spectrum type
        if spectrum_type == "s"
            # Secondary electron spectrum: empirical formula for atomic oxygen
            E_excess = E_primary_energy - E_ionization_threshold
            secondary_spectrum = B_factor ./ (1 .+ (E_grid ./ A_factor).^(5/3)) .* (E_grid .< E_excess / 2)
            return secondary_spectrum

        elseif spectrum_type == "c"
            # Cascading (degraded primary) spectrum: extract directly from transfer matrix
            i_threshold = findmin(abs.(E_ionization_threshold .- ionization_thresholds))[2]
            i_primary = findmin(abs.(E_grid_for_Q .- E_primary_energy))[2]
            cascading_spectrum = Q_transfer_matrix[i_primary, 1:length(E_grid), i_threshold]
            return cascading_spectrum
        end
    end
end










# ======================================================================================== #
#                                CALCULATION FUNCTIONS                                    #
# ======================================================================================== #

# Notes:
# Each function has a version where the integration is done with hcubature, and one
# version where it is done with quadgk. The integration might be easier to understand when
# looking at the quadgk versions.
# From my testing, the quadgk versions are faster for N2 and O2, but slower for O (I believe
# because of the ^(5/3) factor that get called more times in the quadgk version). The
# hcubature versions are the one used by default because they were the one used before the
# refactor of the cascading functions.

"""
    calculate_cascading_N2(E_grid, dE, lorentzian_width)

Calculate cascading transfer matrices for N₂ ionization.

# Arguments
- `E_grid`: Energy grid for electrons (eV)
- `dE`: Energy grid step size (eV)
- `lorentzian_width`: Width of the Lorentzian distribution (eV)

# Returns
- `Q_transfer_matrix`: Transfer matrix [n_energies, n_energies, n_thresholds]

# Reference
Equation from Itikawa 1986 J. Phys. Chem. Ref. Data
"""
function calculate_cascading_N2(E_grid, dE, lorentzian_width = 11.4)
    # N₂ ionization thresholds for different states (eV)
    ionization_thresholds = [15.581, 16.73, 18.75, 24, 42]

    n_energies = length(E_grid)
    n_thresholds = length(ionization_thresholds)
    Q_transfer_matrix = zeros(n_energies, n_energies, n_thresholds)

    println("Pre-calculating all energy-degradations for e⁻ - N₂ ionizing collisions.")
    println("Starting at ", Dates.format(now(), "HH:MM:SS"))

    # Loop over ionization thresholds (reverse order for display purposes)
    for i_threshold in n_thresholds:-1:1
        threshold = ionization_thresholds[i_threshold]

        # Find first energy bin above ionization threshold
        i_min_primary = findfirst(x -> x > threshold, E_grid)

        progress = Progress(n_energies - i_min_primary + 1;
                           desc=string("N₂ threshold ", n_thresholds - i_threshold + 1, "/", n_thresholds),
                           color=:blue, dt=0.01)

        # Loop over primary electron energy bins
        for i_primary in i_min_primary:n_energies
            E_primary_bin_min = E_grid[i_primary]
            E_primary_bin_max = E_grid[i_primary] + dE[i_primary]

            # Maximum degraded energy: half of excess energy
            E_degraded_max_physical = (E_primary_bin_min - threshold) / 2
            i_max_degraded = findlast(x -> x < E_degraded_max_physical, E_grid)

            if !isnothing(i_max_degraded)
                # Loop over degraded electron energy bins
                for i_degraded in i_max_degraded:(i_primary - 1)
                    E_degraded_bin_min = E_grid[i_degraded]
                    E_degraded_bin_max = min(E_grid[i_degraded] + dE[i_degraded],
                                             E_primary_bin_max - threshold)

                    # Integrate only if limits are physical
                    if E_degraded_bin_max > E_degraded_bin_min
                        function integrand(vars)
                            E_degraded, u_primary = vars

                            # To be able to integrate over a cube, we perform a change of
                            # variable. The variable E_primary that varies from E_prim_min
                            # to E_primary_bin_max is replaced by a variable u_primary
                            # that varies from [0, 1].
                            # We have to do this because the lower limit E_prim_min varies with respect
                            # to the other variable of our integral (E-degraded). As such, its
                            # integration domain has the shape of a triangle/trapezoid.
                            E_prim_min = max(E_primary_bin_min, E_degraded + threshold)
                            E_primary_mapped = E_prim_min + u_primary * (E_primary_bin_max - E_prim_min)
                            # We need a jacobian due to our change of variable
                            jacobian = E_primary_bin_max - E_prim_min

                            # Calculate secondary electron energy
                            E_secondary = E_primary_mapped - threshold - E_degraded
                            # Lorentzian distribution
                            lorentzian = 1.0 / (E_secondary^2 + lorentzian_width^2)

                            return lorentzian * jacobian
                        end

                        result, _ = hcubature(integrand,
                                             [E_degraded_bin_min, 0.0],
                                             [E_degraded_bin_max, 1.0],
                                             rtol=1e-4, atol=1e-10)
                        Q_transfer_matrix[i_primary, i_degraded, i_threshold] = result
                    end
                end
            end
            next!(progress)
        end
    end

    return Q_transfer_matrix, E_grid, ionization_thresholds
end


"""
    calculate_cascading_N2_quadgk(E_grid, dE, lorentzian_width)

Calculate cascading transfer matrices for N₂ ionization using nested quadgk integrals.

This is an alternative implementation to `calculate_cascading_N2` that uses two nested
one-dimensional integrals instead of a hypercube transformation. The integration is performed
over the physical domain directly:
- Outer integral: over degraded electron energy E_degraded
- Inner integral: over primary electron energy E_primary

# Arguments
- `E_grid`: Energy grid for electrons (eV)
- `dE`: Energy grid step size (eV)
- `lorentzian_width`: Width of the Lorentzian distribution (eV)

# Returns
- `Q_transfer_matrix`: Transfer matrix [n_energies, n_energies, n_thresholds]
- `E_grid`: Energy grid (returned for consistency)
- `ionization_thresholds`: Array of ionization threshold energies

# Reference
Equation from Itikawa 1986 J. Phys. Chem. Ref. Data
"""
function calculate_cascading_N2_quadgk(E_grid, dE, lorentzian_width = 11.4)
    # N₂ ionization thresholds for different states (eV)
    ionization_thresholds = [15.581, 16.73, 18.75, 24, 42]

    n_energies = length(E_grid)
    n_thresholds = length(ionization_thresholds)
    Q_transfer_matrix = zeros(n_energies, n_energies, n_thresholds)

    println("Pre-calculating all energy-degradations for e⁻ - N₂ ionizing collisions (quadgk version).")
    println("Starting at ", Dates.format(now(), "HH:MM:SS"))

    # Loop over ionization thresholds (reverse order for display purposes)
    for i_threshold in n_thresholds:-1:1
        threshold = ionization_thresholds[i_threshold]

        # Find first energy bin above ionization threshold
        i_min_primary = findfirst(x -> x > threshold, E_grid)

        progress = Progress(n_energies - i_min_primary + 1;
                           desc=string("N₂ threshold (quadgk) ", n_thresholds - i_threshold + 1, "/", n_thresholds),
                           color=:blue, dt=0.01)

        # Loop over primary electron energy bins
        Threads.@threads for i_primary in i_min_primary:n_energies
            E_primary_bin_min = E_grid[i_primary]
            E_primary_bin_max = E_grid[i_primary] + dE[i_primary]

            # Maximum degraded energy: half of excess energy
            E_degraded_max_physical = (E_primary_bin_min - threshold) / 2
            i_max_degraded = findlast(x -> x < E_degraded_max_physical, E_grid)

            if !isnothing(i_max_degraded)
                # Loop over degraded electron energy bins
                for i_degraded in i_max_degraded:(i_primary - 1)
                    E_degraded_bin_min = E_grid[i_degraded]
                    E_degraded_bin_max = min(E_grid[i_degraded] + dE[i_degraded],
                                             E_primary_bin_max - threshold)
                    if E_degraded_bin_max > E_degraded_bin_min
                        # Outer integral over degraded electron energy
                        result, _ = quadgk(E_degraded_bin_min, E_degraded_bin_max) do E_degraded
                            # For each E_degraded, determine the valid range of E_primary
                            E_primary_min_local = max(E_primary_bin_min, E_degraded + threshold)
                            E_primary_max_local = E_primary_bin_max

                            # Inner integral over primary electron energy
                            if E_primary_max_local > E_primary_min_local
                                inner_result, _ = quadgk(E_primary_min_local, E_primary_max_local) do E_primary
                                    # Calculate secondary electron energy
                                    E_secondary = E_primary - threshold - E_degraded
                                    # Lorentzian distribution for secondary electron
                                    lorentzian = 1.0 / (E_secondary^2 + lorentzian_width^2)
                                    return lorentzian
                                end
                                return inner_result
                            else
                                return 0.0
                            end
                        end

                        Q_transfer_matrix[i_primary, i_degraded, i_threshold] = result
                    end
                end
            end
            next!(progress)
        end
    end

    return Q_transfer_matrix, E_grid, ionization_thresholds
end

"""
    calculate_cascading_O2(E_grid, dE, lorentzian_width)

Calculate cascading transfer matrices for O₂ ionization.

# Arguments
- `E_grid`: Energy grid for electrons (eV)
- `dE`: Energy grid step size (eV)
- `lorentzian_width`: Width of the Lorentzian distribution (eV)

# Returns
- `Q_transfer_matrix`: Transfer matrix [n_energies, n_energies, n_thresholds]
"""
function calculate_cascading_O2(E_grid, dE, lorentzian_width = 15.2)
    # O₂ ionization thresholds (eV)
    ionization_thresholds = [12.072, 16.1, 16.9, 18.2, 18.9, 32.51]

    n_energies = length(E_grid)
    n_thresholds = length(ionization_thresholds)
    Q_transfer_matrix = zeros(n_energies, n_energies, n_thresholds)

    println("Pre-calculating all energy-degradations for e⁻ - O₂ ionizing collisions.")
    println("Starting at ", Dates.format(now(), "HH:MM:SS"))

    for i_threshold in n_thresholds:-1:1
        threshold = ionization_thresholds[i_threshold]
        i_min_primary = findfirst(x -> x > threshold, E_grid)

        progress = Progress(n_energies - i_min_primary + 1;
                           desc=string("O₂ threshold ", n_thresholds - i_threshold + 1, "/", n_thresholds),
                           color=:blue, dt=0.01)

        for i_primary in i_min_primary:n_energies
            E_primary_bin_min = E_grid[i_primary]
            E_primary_bin_max = E_grid[i_primary] + dE[i_primary]

            E_degraded_max_physical = (E_primary_bin_min - threshold) / 2
            i_max_degraded = findlast(x -> x < E_degraded_max_physical, E_grid)

            if !isnothing(i_max_degraded)
                for i_degraded in i_max_degraded:(i_primary - 1)
                    E_degraded_bin_min = E_grid[i_degraded]
                    E_degraded_bin_max = min(E_grid[i_degraded] + dE[i_degraded],
                                             E_primary_bin_max - threshold)

                    # Integrate only if limits are physical
                    if E_degraded_bin_max > E_degraded_bin_min
                        function integrand(vars)
                            E_degraded, u_primary = vars

                            # To be able to integrate over a cube, we perform a change of
                            # variable. The variable E_primary that varies from E_prim_min
                            # to E_primary_bin_max is replaced by a variable u_primary
                            # that varies from [0, 1].
                            # We have to do this because the lower limit E_prim_min varies with respect
                            # to the other variable of our integral (E-degraded). As such, its
                            # integration domain has the shape of a triangle/trapezoid.
                            E_prim_min = max(E_primary_bin_min, E_degraded + threshold)
                            E_primary_mapped = E_prim_min + u_primary * (E_primary_bin_max - E_prim_min)
                            # We need a jacobian due to our change of variable
                            jacobian = E_primary_bin_max - E_prim_min

                            # Calculate secondary electron energy
                            E_secondary = E_primary_mapped - threshold - E_degraded
                            # Lorentzian distribution
                            lorentzian = 1.0 / (E_secondary^2 + lorentzian_width^2)

                            return lorentzian * jacobian
                        end

                        result, _ = hcubature(integrand,
                                             [E_degraded_bin_min, 0.0],
                                             [E_degraded_bin_max, 1.0],
                                             rtol=1e-4, atol = 1e-10)
                        Q_transfer_matrix[i_primary, i_degraded, i_threshold] = result
                    end
                end
            end
            next!(progress)
        end
    end

    return Q_transfer_matrix, E_grid, ionization_thresholds
end


"""
    calculate_cascading_O2_quadgk(E_grid, dE, lorentzian_width)

Calculate cascading transfer matrices for O₂ ionization using nested quadgk integrals.

This is an alternative implementation to `calculate_cascading_O2` that uses two nested
one-dimensional integrals instead of a hypercube transformation. The integration is performed
over the physical domain directly:
- Outer integral: over degraded electron energy E_degraded
- Inner integral: over primary electron energy E_primary

# Arguments
- `E_grid`: Energy grid for electrons (eV)
- `dE`: Energy grid step size (eV)
- `lorentzian_width`: Width of the Lorentzian distribution (eV)

# Returns
- `Q_transfer_matrix`: Transfer matrix [n_energies, n_energies, n_thresholds]
- `E_grid`: Energy grid (returned for consistency)
- `ionization_thresholds`: Array of ionization threshold energies
"""
function calculate_cascading_O2_quadgk(E_grid, dE, lorentzian_width = 15.2)
    # O₂ ionization thresholds (eV)
    ionization_thresholds = [12.072, 16.1, 16.9, 18.2, 18.9, 32.51]

    n_energies = length(E_grid)
    n_thresholds = length(ionization_thresholds)
    Q_transfer_matrix = zeros(n_energies, n_energies, n_thresholds)

    println("Pre-calculating all energy-degradations for e⁻ - O₂ ionizing collisions (quadgk version).")
    println("Starting at ", Dates.format(now(), "HH:MM:SS"))

    # Loop over ionization thresholds (reverse order for display purposes)
    for i_threshold in n_thresholds:-1:1
        threshold = ionization_thresholds[i_threshold]

        # Find first energy bin above ionization threshold
        i_min_primary = findfirst(x -> x > threshold, E_grid)

        progress = Progress(n_energies - i_min_primary + 1;
                           desc=string("O₂ threshold (quadgk) ", n_thresholds - i_threshold + 1, "/", n_thresholds),
                           color=:blue, dt=0.01)

        # Loop over primary electron energy bins
        for i_primary in i_min_primary:n_energies
            E_primary_bin_min = E_grid[i_primary]
            E_primary_bin_max = E_grid[i_primary] + dE[i_primary]

            # Maximum degraded energy: half of excess energy
            E_degraded_max_physical = (E_primary_bin_min - threshold) / 2
            i_max_degraded = findlast(x -> x < E_degraded_max_physical, E_grid)

            if !isnothing(i_max_degraded)
                # Loop over degraded electron energy bins
                for i_degraded in i_max_degraded:(i_primary - 1)
                    E_degraded_bin_min = E_grid[i_degraded]
                    E_degraded_bin_max = min(E_grid[i_degraded] + dE[i_degraded],
                                             E_primary_bin_max - threshold)
                    if E_degraded_bin_max > E_degraded_bin_min
                        # Outer integral over degraded electron energy
                        result, _ = quadgk(E_degraded_bin_min, E_degraded_bin_max) do E_degraded
                            # For each E_degraded, determine the valid range of E_primary
                            E_primary_min_local = max(E_primary_bin_min, E_degraded + threshold)
                            E_primary_max_local = E_primary_bin_max

                            # Inner integral over primary electron energy
                            if E_primary_max_local > E_primary_min_local
                                inner_result, _ = quadgk(E_primary_min_local, E_primary_max_local) do E_primary
                                    # Calculate secondary electron energy
                                    E_secondary = E_primary - threshold - E_degraded
                                    # Lorentzian distribution for secondary electron
                                    lorentzian = 1.0 / (E_secondary^2 + lorentzian_width^2)
                                    return lorentzian
                                end
                                return inner_result
                            else
                                return 0.0
                            end
                        end

                        Q_transfer_matrix[i_primary, i_degraded, i_threshold] = result
                    end
                end
            end
            next!(progress)
        end
    end

    return Q_transfer_matrix, E_grid, ionization_thresholds
end


"""
    interpolate_O_parameters(E_primary)

Interpolate energy-dependent parameters for atomic O ionization.

# Arguments
- `E_primary`: Primary electron energy (eV)

# Returns
- `(A_factor, B_factor)`: Tuple of interpolated parameters
"""
function interpolate_O_parameters(E_primary)
    energy_params = [100, 200, 500, 1000, 2000]  # eV
    B_params = [7.18, 4.97, 2.75, 1.69, 1.02] .* 1e-22
    A_params = [12.6, 13.7, 14.1, 14.0, 13.7]

    if (minimum(energy_params) < E_primary < maximum(energy_params))
        A_factor = linear_interpolation(energy_params, A_params)(E_primary)
        B_factor = linear_interpolation(energy_params, B_params)(E_primary)
    elseif E_primary <= minimum(energy_params)
        A_factor = A_params[1]
        B_factor = B_params[1]
    else
        A_factor = A_params[end]
        B_factor = B_params[end]
    end

    return (A_factor, B_factor)
end


"""
    calculate_cascading_O(E_grid, dE)

Calculate cascading transfer matrices for atomic O ionization.

# Arguments
- `E_grid`: Energy grid for electrons (eV)
- `dE`: Energy grid step size (eV)

# Returns
- `Q_transfer_matrix`: Transfer matrix [n_energies, n_energies, n_thresholds]
"""
function calculate_cascading_O(E_grid, dE)
    # O ionization thresholds for different states (eV)
    ionization_thresholds = [13.618, 16.9, 18.6, 28.5]

    n_energies = length(E_grid)
    n_thresholds = length(ionization_thresholds)
    Q_transfer_matrix = zeros(n_energies, n_energies, n_thresholds)

    println("Pre-calculating all energy-degradations for e⁻ - O ionizing collisions.")
    println("Starting at ", Dates.format(now(), "HH:MM:SS"))

    # Loop over ionization thresholds (reverse order for display purposes)
    for i_threshold in n_thresholds:-1:1
        threshold = ionization_thresholds[i_threshold]

        # Find first energy bin above ionization threshold
        i_min_primary = findfirst(x -> x > threshold, E_grid)

        progress = Progress(n_energies - i_min_primary + 1;
                           desc=string("O threshold ", n_thresholds - i_threshold + 1, "/", n_thresholds),
                           color=:blue, dt=0.01)

        # Loop over primary electron energy bins
        for i_primary in i_min_primary:n_energies
            # Update energy-dependent parameters for this primary energy bin
            E_primary_bin = E_grid[i_primary]
            A_local, B_local = interpolate_O_parameters(E_primary_bin)

            E_primary_bin_min = E_grid[i_primary]
            E_primary_bin_max = E_grid[i_primary] + dE[i_primary]

            # Maximum degraded energy: half of excess energy
            E_degraded_max_physical = (E_primary_bin_min - threshold) / 2
            i_max_degraded = findlast(x -> x < E_degraded_max_physical, E_grid)

            if !isnothing(i_max_degraded)
                # Loop over degraded electron energy bins
                for i_degraded in i_max_degraded:(i_primary - 1)
                    E_degraded_bin_min = E_grid[i_degraded]
                    E_degraded_bin_max = min(E_grid[i_degraded] + dE[i_degraded],
                                             E_primary_bin_max - threshold)

                    # Integrate only if limits are physical
                    if E_degraded_bin_max > E_degraded_bin_min
                        function integrand(vars)
                            E_degraded, u_primary = vars

                            # To be able to integrate over a cube, we perform a change of
                            # variable. The variable E_primary that varies from E_prim_min
                            # to E_primary_bin_max is replaced by a variable u_primary
                            # that varies from [0, 1].
                            # We have to do this because the lower limit E_prim_min varies with respect
                            # to the other variable of our integral (E-degraded). As such, its
                            # integration domain has the shape of a triangle/trapezoid.
                            E_prim_min = max(E_primary_bin_min, E_degraded + threshold)
                            E_primary_mapped = E_prim_min + u_primary * (E_primary_bin_max - E_prim_min)
                            # We need a jacobian due to our change of variable
                            jacobian = E_primary_bin_max - E_prim_min

                            # Calculate secondary electron energy
                            E_secondary = E_primary_mapped - threshold - E_degraded
                            # Empirical distribution for atomic oxygen
                            distribution = B_local / (1 + (E_secondary / A_local)^(5/3))

                            return distribution * jacobian
                        end

                        result, _ = hcubature(integrand,
                                             [E_degraded_bin_min, 0.0],
                                             [E_degraded_bin_max, 1.0],
                                             rtol=1e-4, atol=0)
                        Q_transfer_matrix[i_primary, i_degraded, i_threshold] = result
                    end
                end
            end
            next!(progress)
        end
    end

    return Q_transfer_matrix, E_grid, ionization_thresholds
end


"""
    calculate_cascading_O_quadgk(E_grid, dE)

Calculate cascading transfer matrices for atomic O ionization using nested quadgk integrals.

This is an alternative implementation to `calculate_cascading_O` that uses two nested
one-dimensional integrals instead of a hypercube transformation. The integration is performed
over the physical domain directly:
- Outer integral: over degraded electron energy E_degraded
- Inner integral: over primary electron energy E_primary

# Arguments
- `E_grid`: Energy grid for electrons (eV)
- `dE`: Energy grid step size (eV)

# Returns
- `Q_transfer_matrix`: Transfer matrix [n_energies, n_energies, n_thresholds]
- `E_grid`: Energy grid (returned for consistency)
- `ionization_thresholds`: Array of ionization threshold energies
"""
function calculate_cascading_O_quadgk(E_grid, dE)
    # O ionization thresholds for different states (eV)
    ionization_thresholds = [13.618, 16.9, 18.6, 28.5]

    n_energies = length(E_grid)
    n_thresholds = length(ionization_thresholds)
    Q_transfer_matrix = zeros(n_energies, n_energies, n_thresholds)

    println("Pre-calculating all energy-degradations for e⁻ - O ionizing collisions (quadgk version).")
    println("Starting at ", Dates.format(now(), "HH:MM:SS"))

    # Loop over ionization thresholds (reverse order for display purposes)
    for i_threshold in n_thresholds:-1:1
        threshold = ionization_thresholds[i_threshold]

        # Find first energy bin above ionization threshold
        i_min_primary = findfirst(x -> x > threshold, E_grid)

        progress = Progress(n_energies - i_min_primary + 1;
                           desc=string("O threshold (quadgk) ", n_thresholds - i_threshold + 1, "/", n_thresholds),
                           color=:blue, dt=0.01)

        # Loop over primary electron energy bins
        for i_primary in i_min_primary:n_energies
            # Update energy-dependent parameters for this primary energy bin
            E_primary_bin = E_grid[i_primary]
            A_local, B_local = interpolate_O_parameters(E_primary_bin)

            E_primary_bin_min = E_grid[i_primary]
            E_primary_bin_max = E_grid[i_primary] + dE[i_primary]

            # Maximum degraded energy: half of excess energy
            E_degraded_max_physical = (E_primary_bin_min - threshold) / 2
            i_max_degraded = findlast(x -> x < E_degraded_max_physical, E_grid)

            if !isnothing(i_max_degraded)
                # Loop over degraded electron energy bins
                for i_degraded in i_max_degraded:(i_primary - 1)
                    E_degraded_bin_min = E_grid[i_degraded]
                    E_degraded_bin_max = min(E_grid[i_degraded] + dE[i_degraded],
                                             E_primary_bin_max - threshold)
                    if E_degraded_bin_max > E_degraded_bin_min
                        # Outer integral over degraded electron energy
                        result, _ = quadgk(E_degraded_bin_min, E_degraded_bin_max) do E_degraded
                            # For each E_degraded, determine the valid range of E_primary
                            E_primary_min_local = max(E_primary_bin_min, E_degraded + threshold)
                            E_primary_max_local = E_primary_bin_max

                            # Inner integral over primary electron energy
                            if E_primary_max_local > E_primary_min_local
                                inner_result, _ = quadgk(E_primary_min_local, E_primary_max_local) do E_primary
                                    # Calculate secondary electron energy
                                    E_secondary = E_primary - threshold - E_degraded
                                    # Empirical distribution for atomic oxygen
                                    distribution = B_local / (1 + (E_secondary / A_local)^(5/3))
                                    return distribution
                                end
                                return inner_result
                            else
                                return 0.0
                            end
                        end

                        Q_transfer_matrix[i_primary, i_degraded, i_threshold] = result
                    end
                end
            end
            next!(progress)
        end
    end

    return Q_transfer_matrix, E_grid, ionization_thresholds
end
