using Dates: Dates, now
using HCubature: hcubature
using Interpolations: linear_interpolation
using MAT: matopen


# ======================================================================================== #
#                                   UTILITY FUNCTIONS                                      #
# ======================================================================================== #


"""
    find_cascading_file(E_edges, species_dir)

Search for a pre-computed cascading spectra file with matching energy grid.

# Arguments
- `E_edges`: Energy grid edges to match
- `species_dir`: Directory containing cascading data files for the species

# Returns
- `(file_found, filepath)`: Tuple of boolean and filepath string
"""
function find_cascading_file(E_edges, species_dir)
    cascading_files = readdir(species_dir)

    for filename in cascading_files
        if !isdir(filename)
            try
                filepath = joinpath(species_dir, filename)
                file = matopen(filepath)
                E_edges_saved = read(file, "E_edges")
                close(file)

                # Check if saved grid matches requested grid
                if length(E_edges) <= length(E_edges_saved) && E_edges_saved[1:length(E_edges)] == E_edges
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
- `(Q_transfer_matrix, E_edges_for_Q, ionization_thresholds)`: Tuple of loaded data
"""
function load_cascading_matrices(filepath)
    println("Loading cascading-matrices from file: ", basename(filepath))
    file = matopen(filepath)
    Q_transfer_matrix = read(file, "Q")
    ionization_thresholds = read(file, "E_ionizations")
    E_edges_for_Q = read(file, "E_edges")
    close(file)

    return (Q_transfer_matrix, E_edges_for_Q, ionization_thresholds)
end


"""
    save_cascading_matrices(Q_transfer_matrix, E_edges_for_Q, ionization_thresholds, species_dir, species_name)

Save calculated cascading matrices to a file.

# Arguments
- `Q_transfer_matrix`: Transfer matrix to save
- `E_edges_for_Q`: Energy grid edges used for calculations
- `ionization_thresholds`: Ionization threshold energies
- `species_dir`: Directory to save the file
- `species_name`: Name of the species (for filename)
"""
function save_cascading_matrices(Q_transfer_matrix, E_edges_for_Q, ionization_thresholds, species_dir, species_name)
    filename = joinpath(species_dir,
                       string("CascadingSpec", species_name, "ionization_",
                             Dates.format(now(), "yyyymmdd-HHMMSS"),
                             ".mat"))
    file = matopen(filename, "w")
    write(file, "Q", Q_transfer_matrix)
    write(file, "E_edges", E_edges_for_Q)
    write(file, "E_ionizations", ionization_thresholds)
    close(file)
end


mutable struct SpeciesCascadingCache
    species::Symbol
    Q_transfer_matrix::Array{Float64, 3}
    E_edges_for_Q::Vector{Float64}
    ionization_thresholds::Vector{Float64}
end

function SpeciesCascadingCache(species::Symbol)
    return SpeciesCascadingCache(species, zeros(0, 0, 0), Float64[], Float64[])
end

struct CascadingCache
    species::NTuple{3, SpeciesCascadingCache}
end

function CascadingCache()
    return CascadingCache((
        SpeciesCascadingCache(:N2),
        SpeciesCascadingCache(:O2),
        SpeciesCascadingCache(:O),
    ))
end

Base.getindex(cache::CascadingCache, index::Int) = cache.species[index]
Base.length(cache::CascadingCache) = length(cache.species)
Base.iterate(cache::CascadingCache, state...) = iterate(cache.species, state...)

species_label(cache::SpeciesCascadingCache) = String(cache.species)

function cascading_lorentzian_width(cache::SpeciesCascadingCache)
    if cache.species == :N2
        return 11.4
    elseif cache.species == :O2
        return 15.2
    end

    error("Lorentzian width is only defined for N2 and O2 cascading.")
end

cascading_species_dir(cache::SpeciesCascadingCache) =
    pkgdir(AURORA, "internal_data", "e_cascading", species_label(cache))

function needs_cascading_reload(cache::SpeciesCascadingCache, energy_grid::EnergyGrid)
    E_edges = energy_grid.E_edges
    return isempty(cache.Q_transfer_matrix) ||
           length(E_edges) > length(cache.E_edges_for_Q) ||
           cache.E_edges_for_Q[1:length(E_edges)] != E_edges
end

function calculate_cascading_matrices(species::Symbol, energy_grid::EnergyGrid)
    if species == :N2
        return calculate_cascading_N2(energy_grid, 11.4)
    elseif species == :O2
        return calculate_cascading_O2(energy_grid, 15.2)
    elseif species == :O
        return calculate_cascading_O(energy_grid)
    end

    error("Unsupported cascading species: $(species)")
end

function ensure_cascading_loaded!(cache::SpeciesCascadingCache, energy_grid::EnergyGrid;
                                  force_recompute::Bool = false)
    force_recompute || needs_cascading_reload(cache, energy_grid) || return cache

    species_dir = cascading_species_dir(cache)
    file_found, filepath = force_recompute ? (false, "") : find_cascading_file(energy_grid.E_edges, species_dir)

    if file_found
        cache.Q_transfer_matrix, cache.E_edges_for_Q, cache.ionization_thresholds =
            load_cascading_matrices(filepath)
    else
        if force_recompute
            println("Forcing recomputation of cascading-matrices (ignoring cached files on disk).")
        else
            println("Could not find a file with matching energy grid.")
        end

        cache.Q_transfer_matrix, cache.E_edges_for_Q, cache.ionization_thresholds =
            calculate_cascading_matrices(cache.species, energy_grid)

        save_cascading_matrices(cache.Q_transfer_matrix, cache.E_edges_for_Q,
                                cache.ionization_thresholds, species_dir,
                                species_label(cache))
    end

    return cache
end

function secondary_spectrum(cache::SpeciesCascadingCache, energy_grid::EnergyGrid,
                            E_primary_energy, E_ionization_threshold)
    ensure_cascading_loaded!(cache, energy_grid)

    E_left = @view(energy_grid.E_edges[1:end-1])
    E_excess = E_primary_energy - E_ionization_threshold

    if cache.species == :O
        A_factor, B_factor = interpolate_O_parameters(E_primary_energy)
        # A_factor *= 1.25
        return B_factor ./ (1 .+ (E_left ./ A_factor).^(5 / 3)) .* (E_left .< E_excess / 2)
    end

    lorentzian_width = cascading_lorentzian_width(cache)
    return 1 ./ (lorentzian_width^2 .+ E_left.^2) .* (E_left .< E_excess / 2)
end

function primary_spectrum(cache::SpeciesCascadingCache, energy_grid::EnergyGrid,
                          E_primary_energy, E_ionization_threshold)
    ensure_cascading_loaded!(cache, energy_grid)

    i_threshold = findmin(abs.(E_ionization_threshold .- cache.ionization_thresholds))[2]
    i_primary = findmin(abs.(@view(cache.E_edges_for_Q[1:end-1]) .- E_primary_energy))[2]
    return cache.Q_transfer_matrix[i_primary, 1:energy_grid.n, i_threshold]
end







# ======================================================================================== #
#                                CALCULATION FUNCTIONS                                     #
# ======================================================================================== #
"""
    calculate_cascading_N2(energy_grid, lorentzian_width)

Calculate cascading transfer matrices for N₂ ionization.

# Arguments
- `energy_grid`: Energy grid
- `lorentzian_width`: Width of the Lorentzian distribution (eV)

# Returns
- `Q_transfer_matrix`: Transfer matrix [n_E, n_E, n_thresholds]

# Reference
Equation from Itikawa 1986 J. Phys. Chem. Ref. Data
"""
function calculate_cascading_N2(energy_grid::EnergyGrid, lorentzian_width = 11.4; verbose = true)
    E_edges = energy_grid.E_edges
    E_left = @view(E_edges[1:end-1])
    n_E = energy_grid.n

    # N₂ ionization thresholds for different states (eV)
    ionization_thresholds = [15.581, 16.73, 18.75, 24, 42]

    n_thresholds = length(ionization_thresholds)
    Q_transfer_matrix = zeros(n_E, n_E, n_thresholds)

    verbose && print("Calculating energy-degradation transfer matrices for e⁻ - N₂ ionizing collisions...")

    # Loop over ionization thresholds (reverse order for display purposes)
    for i_threshold in n_thresholds:-1:1
        threshold = ionization_thresholds[i_threshold]

        # Find first energy bin above ionization threshold
        i_min_primary = findfirst(x -> x > threshold, E_left)

        # Loop over primary electron energy bins
        for i_primary in i_min_primary:n_E
            E_primary_bin_min = E_edges[i_primary]
            E_primary_bin_max = E_edges[i_primary + 1]

            # Maximum degraded energy: half of excess energy
            E_degraded_max_physical = (E_primary_bin_min - threshold) / 2
            i_max_degraded = findlast(x -> x < E_degraded_max_physical, E_left)

            if !isnothing(i_max_degraded)
                # Loop over degraded electron energy bins
                for i_degraded in i_max_degraded:(i_primary - 1)
                    E_degraded_bin_min = E_edges[i_degraded]
                    E_degraded_bin_max = min(E_edges[i_degraded + 1],
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
                                             [E_degraded_bin_max, 1.0])
                        Q_transfer_matrix[i_primary, i_degraded, i_threshold] = result
                    end
                end
            end
        end
    end

    verbose && println(" done.")
    return Q_transfer_matrix, E_edges, ionization_thresholds
end

"""
    calculate_cascading_O2(energy_grid, lorentzian_width)

Calculate cascading transfer matrices for O₂ ionization.

# Arguments
- `energy_grid`: Energy grid
- `lorentzian_width`: Width of the Lorentzian distribution (eV)

# Returns
- `Q_transfer_matrix`: Transfer matrix [n_E, n_E, n_thresholds]
"""
function calculate_cascading_O2(energy_grid::EnergyGrid, lorentzian_width = 15.2; verbose = true)
    E_edges = energy_grid.E_edges
    E_left = @view(E_edges[1:end-1])
    n_E = energy_grid.n

    # O₂ ionization thresholds (eV)
    ionization_thresholds = [12.072, 16.1, 16.9, 18.2, 18.9, 32.51]

    n_thresholds = length(ionization_thresholds)
    Q_transfer_matrix = zeros(n_E, n_E, n_thresholds)

    verbose && print("Calculating energy-degradation transfer matrices for e⁻ - O₂ ionizing collisions...")

    for i_threshold in n_thresholds:-1:1
        threshold = ionization_thresholds[i_threshold]
        i_min_primary = findfirst(x -> x > threshold, E_left)

        for i_primary in i_min_primary:n_E
            E_primary_bin_min = E_edges[i_primary]
            E_primary_bin_max = E_edges[i_primary + 1]

            E_degraded_max_physical = (E_primary_bin_min - threshold) / 2
            i_max_degraded = findlast(x -> x < E_degraded_max_physical, E_left)

            if !isnothing(i_max_degraded)
                for i_degraded in i_max_degraded:(i_primary - 1)
                    E_degraded_bin_min = E_edges[i_degraded]
                    E_degraded_bin_max = min(E_edges[i_degraded + 1],
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
                                             [E_degraded_bin_max, 1.0])
                        Q_transfer_matrix[i_primary, i_degraded, i_threshold] = result
                    end
                end
            end
        end
    end

    verbose && println(" done.")
    return Q_transfer_matrix, E_edges, ionization_thresholds
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
    calculate_cascading_O(energy_grid)

Calculate cascading transfer matrices for atomic O ionization.

# Arguments
- `energy_grid`: Energy grid

# Returns
- `Q_transfer_matrix`: Transfer matrix [n_E, n_E, n_thresholds]
"""
function calculate_cascading_O(energy_grid::EnergyGrid; verbose = true)
    E_edges = energy_grid.E_edges
    E_left = @view(E_edges[1:end-1])
    n_E = energy_grid.n

    # O ionization thresholds for different states (eV)
    ionization_thresholds = [13.618, 16.9, 18.6, 28.5]

    n_thresholds = length(ionization_thresholds)
    Q_transfer_matrix = zeros(n_E, n_E, n_thresholds)

    verbose && print("Calculating energy-degradation transfer matrices for e⁻ - O ionizing collisions...")

    # Loop over ionization thresholds (reverse order for display purposes)
    for i_threshold in n_thresholds:-1:1
        threshold = ionization_thresholds[i_threshold]

        # Find first energy bin above ionization threshold
        i_min_primary = findfirst(x -> x > threshold, E_left)

        # Loop over primary electron energy bins
        for i_primary in i_min_primary:n_E
            # Update energy-dependent parameters for this primary energy bin
            E_primary_bin = E_edges[i_primary]
            A_local, B_local = interpolate_O_parameters(E_primary_bin)

            E_primary_bin_min = E_edges[i_primary]
            E_primary_bin_max = E_edges[i_primary + 1]

            # Maximum degraded energy: half of excess energy
            E_degraded_max_physical = (E_primary_bin_min - threshold) / 2
            i_max_degraded = findlast(x -> x < E_degraded_max_physical, E_left)

            if !isnothing(i_max_degraded)
                # Loop over degraded electron energy bins
                for i_degraded in i_max_degraded:(i_primary - 1)
                    E_degraded_bin_min = E_edges[i_degraded]
                    E_degraded_bin_max = min(E_edges[i_degraded + 1],
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
                                             [E_degraded_bin_max, 1.0])
                        Q_transfer_matrix[i_primary, i_degraded, i_threshold] = result
                    end
                end
            end
        end
    end

    verbose && println(" done.")
    return Q_transfer_matrix, E_edges, ionization_thresholds
end
