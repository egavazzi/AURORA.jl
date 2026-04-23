using Dates: Dates, now
using HCubature: hcubature
using Interpolations: linear_interpolation, Flat
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


# ======================================================================================== #
#           CASCADING SPEC — Species specific cascading physics definition                 #
# ======================================================================================== #


struct CascadingSpec{F}
    name::String
    ionization_thresholds::Vector{Float64}
    secondary_law::F   # callable (E_secondary, E_primary_ref) -> Float64
    cache_dir::String
end

function DefaultCascadingSpecN2()
    ionization_thresholds = [15.581, 16.73, 18.75, 24.0, 42.0]
    law = (E_s, E_p) -> 1.0 / (11.4^2 + E_s^2)
    cache_dir = pkgdir(AURORA, "internal_data", "e_cascading", "N2")
    return CascadingSpec("N2", ionization_thresholds, law, cache_dir)
end

function DefaultCascadingSpecO2()
    ionization_thresholds = [12.072, 16.1, 16.9, 18.2, 18.9, 32.51]
    law = (E_s, E_p) -> 1.0 / (15.2^2 + E_s^2)
    cache_dir = pkgdir(AURORA, "internal_data", "e_cascading", "O2")
    return CascadingSpec("O2", ionization_thresholds, law, cache_dir)
end

function DefaultCascadingSpecO()
    ionization_thresholds = [13.618, 16.9, 18.6, 28.5]
    function interpolate_O_parameters(E_primary)
        energy_params = [100, 200, 500, 1000, 2000]  # eV
        B_params = [7.18, 4.97, 2.75, 1.69, 1.02] .* 1e-22
        A_params = [12.6, 13.7, 14.1, 14.0, 13.7]
        A_factor = linear_interpolation(energy_params, A_params, extrapolation_bc = Flat())(E_primary)
        B_factor = linear_interpolation(energy_params, B_params, extrapolation_bc = Flat())(E_primary)
        return (A_factor, B_factor)
    end
    law = function (E_s, E_p)
        A_factor, B_factor = interpolate_O_parameters(E_p)
        return B_factor / (1 + (E_s / A_factor)^(5 / 3))
    end
    cache_dir = pkgdir(AURORA, "internal_data", "e_cascading", "O")
    return CascadingSpec("O", ionization_thresholds, law, cache_dir)
end


# ======================================================================================== #
#                      CASCADING CACHE — Per-species in-memory cache                       #
# ======================================================================================== #

# Species specific cascading cache container
mutable struct SpeciesCascadingCache
    spec::CascadingSpec
    Q_transfer_matrix::Array{Float64, 3}
    E_edges_for_Q::Vector{Float64}
    ionization_thresholds::Vector{Float64}
end

# Intialization constructor
function SpeciesCascadingCache(spec::CascadingSpec)
    return SpeciesCascadingCache(spec, zeros(0, 0, 0), Float64[], Float64[])
end

# Somewhat temporary container to hold all our three species specific caches
# TODO: To be removed when we have fully moved towards a fully modular setup with things
# attached to species inside the model.
struct CascadingCache
    species::NTuple{3, SpeciesCascadingCache}
end
# And its initialization constructor
function CascadingCache()
    return CascadingCache((
        SpeciesCascadingCache(DefaultCascadingSpecN2()),
        SpeciesCascadingCache(DefaultCascadingSpecO2()),
        SpeciesCascadingCache(DefaultCascadingSpecO()),
    ))
end

Base.getindex(cache::CascadingCache, index::Int) = cache.species[index]
Base.length(cache::CascadingCache) = length(cache.species)
Base.iterate(cache::CascadingCache, state...) = iterate(cache.species, state...)

species_label(cache::SpeciesCascadingCache) = cache.spec.name

function needs_cascading_reload(cache::SpeciesCascadingCache, energy_grid::EnergyGrid)
    E_edges = energy_grid.E_edges
    return isempty(cache.Q_transfer_matrix) ||
           length(E_edges) > length(cache.E_edges_for_Q) ||
           cache.E_edges_for_Q[1:length(E_edges)] != E_edges
end

function ensure_cascading_loaded!(cache::SpeciesCascadingCache, energy_grid::EnergyGrid;
                                  force_recompute::Bool = false)
    force_recompute || needs_cascading_reload(cache, energy_grid) || return cache

    file_found, filepath = force_recompute ? (false, "") :
                           find_cascading_file(energy_grid.E_edges, cache.spec.cache_dir)

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
            calculate_transfer_matrix(cache.spec, energy_grid)

        save_cascading_matrices(cache.Q_transfer_matrix, cache.E_edges_for_Q,
                                cache.ionization_thresholds, cache.spec.cache_dir,
                                cache.spec.name)
    end

    return cache
end

function secondary_spectrum(cache::SpeciesCascadingCache, energy_grid::EnergyGrid,
                            E_primary_energy, E_ionization_threshold)
    ensure_cascading_loaded!(cache, energy_grid)

    E_left = @view(energy_grid.E_edges[1:end-1])
    E_excess = E_primary_energy - E_ionization_threshold

    return cache.spec.secondary_law.(E_left, E_primary_energy) .* (E_left .< E_excess / 2)
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
    calculate_transfer_matrix(spec, energy_grid)

Calculate the energy-degradation transfer matrix for a species defined by its `CascadingSpec`.

The outer loop structure is identical for all species. The only species-specific ingredient is
`spec.secondary_law(E_secondary, E_primary_bin_min) -> Float64`, which describes how secondary
electrons distribute in energy given a primary electron at the left edge of the current bin.

The second argument to `secondary_law` is the **bin left edge energy**, not the integrated
`E_primary_mapped`, so that energy-dependent parameters (e.g. A, B for atomic O) are evaluated
once per primary bin — consistent with how the O case was treated historically.

# Arguments
- `spec`: `CascadingSpec` instance with ionization thresholds, secondary distribution law, and cache directory
- `energy_grid`: Energy grid

# Returns
- `(Q_transfer_matrix, E_edges, ionization_thresholds)`: Transfer matrix [n_E, n_E, n_thresholds],
  energy grid edges, and ionization thresholds
"""
function calculate_transfer_matrix(spec::CascadingSpec, energy_grid::EnergyGrid; verbose = true)
    E_edges = energy_grid.E_edges
    E_left = @view(E_edges[1:end-1])
    n_E = energy_grid.n

    ionization_thresholds = spec.ionization_thresholds
    n_thresholds = length(ionization_thresholds)
    Q_transfer_matrix = zeros(n_E, n_E, n_thresholds)

    verbose && print("Calculating energy-degradation transfer matrices for e⁻ - $(spec.name) ionizing collisions...")

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
                            distribution = spec.secondary_law(E_secondary, E_primary_bin_min)

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
