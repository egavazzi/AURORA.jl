using Dates: Dates, now
using HCubature: hcubature, hcubature_buffer
using Interpolations: linear_interpolation, Flat
using MAT: matopen

# ======================================================================================== #
#           CASCADING SPEC — Species specific cascading physics definition                 #
# ======================================================================================== #


struct CascadingSpec{F}
    name::String
    ionization_thresholds::Vector{Float64}
    secondary_law::F   # callable (E_secondary, E_primary) -> Float64
end

function DefaultCascadingSpecN2()
    ionization_thresholds = [15.581, 16.73, 18.75, 24.0, 42.0]
    law = (E_s, E_p) -> 1.0 / (11.4^2 + E_s^2)
    return CascadingSpec("N2", ionization_thresholds, law)
end

function DefaultCascadingSpecO2()
    ionization_thresholds = [12.072, 16.1, 16.9, 18.2, 18.9, 32.51]
    law = (E_s, E_p) -> 1.0 / (15.2^2 + E_s^2)
    return CascadingSpec("O2", ionization_thresholds, law)
end

function DefaultCascadingSpecO()
    ionization_thresholds = [13.618, 16.9, 18.6, 28.5]
    energy_params = [100, 200, 500, 1000, 2000]  # eV
    B_params = [7.18, 4.97, 2.75, 1.69, 1.02] .* 1e-22
    A_params = [12.6, 13.7, 14.1, 14.0, 13.7]
    A_interp = linear_interpolation(energy_params, A_params, extrapolation_bc = Flat())
    B_interp = linear_interpolation(energy_params, B_params, extrapolation_bc = Flat())

    interpolate_O_parameters = E_primary -> (A_interp(E_primary), B_interp(E_primary))

    law = function (E_s, E_p)
        A_factor, B_factor = interpolate_O_parameters(E_p)
        A_factor *= 1.25  # Empirical correction factor
        return B_factor / (1 + (E_s / A_factor)^(5 / 3))
    end
    return CascadingSpec("O", ionization_thresholds, law)
end


# ======================================================================================== #
#                      CASCADING CACHE — Per-species in-memory cache                       #
# ======================================================================================== #

# Species specific cascading cache container
mutable struct SpeciesCascadingCache{S<:CascadingSpec}
    spec::S
    primary_transfer_matrix::Array{Float64, 3}
    secondary_transfer_matrix::Array{Float64, 3}
    E_edges::Vector{Float64}
    ionization_thresholds::Vector{Float64}
end

# Initialization constructor
function SpeciesCascadingCache(spec::S) where {S<:CascadingSpec}
    return SpeciesCascadingCache{S}(spec, zeros(0, 0, 0), zeros(0, 0, 0), Float64[], Float64[])
end

# Somewhat temporary container to hold all our three species specific caches
# TODO: To be removed when we have fully moved towards a fully modular setup with things
# attached to species inside the model.
struct CascadingCache{T<:Tuple}
    species::T
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

function load_or_compute_cascading!(cache::CascadingCache, energy_grid::EnergyGrid;
                                    force_recompute::Bool = false)
    for species_cache in cache
        load_or_compute_cascading!(species_cache, energy_grid; force_recompute)
    end
    return nothing
end

function load_or_compute_cascading!(cache::SpeciesCascadingCache, energy_grid::EnergyGrid;
                                    force_recompute::Bool = false)
    E_edges = energy_grid.E_edges

    file_found, filepath = force_recompute ? (false, "") :
                           find_cascading_file(cache.spec, E_edges)

    if file_found
        try
            cascading_data = load_cascading_matrices(filepath)
            # Trim them to fit the simulation energy grid
            cache.primary_transfer_matrix = cascading_data[1][1:length(E_edges)-1, 1:length(E_edges)-1, :]
            cache.secondary_transfer_matrix = cascading_data[2][1:length(E_edges)-1, 1:length(E_edges)-1, :]
            cache.E_edges = cascading_data[3][1:length(E_edges)]
            cache.ionization_thresholds = cascading_data[4]
            return nothing
        catch
            println("Found incompatible cascading cache file $(basename(filepath)); recomputing.")
        end
    end

    if force_recompute
        println("Forcing recomputation of cascading-matrices (ignoring cached files on disk).")
    else
        if !file_found
            println("Could not find a file with matching energy grid.")
        end
    end

    cascading_data = calculate_cascading_matrices(cache.spec, E_edges)
    cache.primary_transfer_matrix = cascading_data[1]
    cache.secondary_transfer_matrix = cascading_data[2]
    cache.E_edges = cascading_data[3]
    cache.ionization_thresholds = cascading_data[4]

    save_cascading_matrices(cache.primary_transfer_matrix, cache.secondary_transfer_matrix,
                            cache.E_edges, cache.ionization_thresholds,
                            cache.spec.name)

    return nothing
end







# ======================================================================================== #
#                                CALCULATION FUNCTIONS                                     #
# ======================================================================================== #

#=
The implementation here can be a bit hard to follow, so here are some explanations.

For each primary energy bin (outer loop), we want to calculate the distribution of
probabilities of degraded primary. To do this, we loop over allowed degraded primary bins
(inner loop). For each degraded primary bin, we fix a degraded primary energy inside of it,
and integrate over the allowed primary energies that can produce this degraded primary energy.
For each of these allowed primary energies and fixed degraded primary energy, the secondary
energy is determined by energy conservation, and we can evaluate the secondary distribution
law for that specific secondary energy to get the probability of this specific combination
of primary and degraded primary energies. Then we fix another degraded primary energy inside
the degraded bin, and repeat until we have covered the whole degraded primary bin.
Mathematically, all of this translates into a double integral.

We do the same thing for the secondary transfer matrix, except that we now inner loop over
allowed secondary bins, and for each fixed secondary energy we integrate over the allowed
primary energies that can produce this secondary energy.
=#

"""
    DegradedCascadingIntegrand(E_primary_bin_min, E_primary_bin_max, threshold, secondary_law)

Callable integrand used by `hcubature` to integrate the degraded-primary transfer matrix.
"""
struct DegradedCascadingIntegrand{FT, F}
    E_primary_bin_min::FT
    E_primary_bin_max::FT
    threshold::FT
    secondary_law::F
end

"""
    DegradedCascadingIntegrand(E_primary_bin_min, E_primary_bin_max, threshold, secondary_law)

Callable integrand used by `hcubature` to integrate the secondary transfer matrix.
"""
struct SecondaryCascadingIntegrand{FT, F}
    E_primary_bin_min::FT
    E_primary_bin_max::FT
    threshold::FT
    secondary_law::F
end



function (integrand::DegradedCascadingIntegrand)(vars)
    # integration variables passed by hcubature; u_primary ∈ [0, 1] is the mapped primary-energy coordinate
    E_degraded, u_primary = vars

    # The lower bound of E_primary depends on E_degraded (it must exceed E_degraded + threshold
    # so that the secondary electron energy is positive).
    # The upper bound of E_primary also depends on E_degraded (it cannot exceed
    # 2 * E_degraded + threshold, otherwise the secondary would be more energetic than
    # the degraded primary.
    # This comes from
    #       E_secondary = E_primary - threshold - E_degraded
    # and
    #       E_secondary ≤ E_degraded
    #   =>  E_primary - threshold - E_degraded ≤ E_degraded
    #   =>  E_primary ≤ 2 * E_degraded + threshold
    #
    # This makes the 2D integration domain non-rectangular (trapezoidal).
    # To handle this with hcubature, we map E_primary onto u_primary ∈ [0, 1] where
    # u_primary = 0 corresponds to E_primary_lower and u_primary = 1 corresponds to
    # E_primary_upper. The Jacobian dE_primary/du_primary of this transformation is
    # included below.
    E_primary_lower = max(integrand.E_primary_bin_min, E_degraded + integrand.threshold)
    E_primary_upper = min(integrand.E_primary_bin_max, integrand.threshold + 2 * E_degraded)
    jacobian = E_primary_upper - E_primary_lower
    jacobian <= 0 && return 0.0 # not really necessary as
                                # (E_primary_lower ≤ E_primary_upper)
                                # but one is never too safe

    E_primary = E_primary_lower + u_primary * jacobian
    E_secondary = E_primary - integrand.threshold - E_degraded

    return jacobian * integrand.secondary_law(E_secondary, E_primary)
end

function (integrand::SecondaryCascadingIntegrand)(vars)
    E_secondary, u_primary = vars

    # For a fixed secondary energy, we enforce E_secondary <= E_degraded with
    # E_degraded = E_primary - threshold - E_secondary, which implies
    # E_primary >= threshold + 2 * E_secondary.
    E_primary_lower = max(integrand.E_primary_bin_min,
                          integrand.threshold + 2 * E_secondary)
    E_primary_upper = integrand.E_primary_bin_max
    jacobian = E_primary_upper - E_primary_lower
    jacobian <= 0 && return 0.0

    E_primary = E_primary_lower + u_primary * jacobian

    return jacobian * integrand.secondary_law(E_secondary, E_primary)
end


"""
    calculate_cascading_matrices(spec, E_edges; verbose=true)

Calculate the energy-degradation transfer matrix for a species defined by its `CascadingSpec`.
The matrices can later be used to directly get the degraded primary and secondary electron
distributions for any given primary energy and ionization threshold.

The outer loop structure is identical for all species. The only species-specific ingredients are
- `spec.secondary_law(E_secondary, E_primary) -> Float64`, which describes how secondary
    electrons distribute in energy given a primary electron at the current integration energy.
- `spec.ionization_thresholds`, which define the ionization thresholds for the species and thus
    the number of transfer matrices to calculate.

# Arguments
- `spec::CascadingSpec` contains species name, ionization thresholds and secondary distribution law
- `E_edges`: Energy grid edges to match (eV)

# Returns
- `(primary_transfer_matrix, secondary_transfer_matrix, E_edges, ionization_thresholds)`:
    degraded-primary transfer matrix [n_E, n_E, n_thresholds], secondary transfer matrix
    [n_E, n_E, n_thresholds], energy grid edges, and ionization thresholds
"""
function calculate_cascading_matrices(spec::CascadingSpec, E_edges; verbose = true)
    E_left = @view(E_edges[1:end-1])
    n_E = length(E_left) # number of energy bins is one less than number of edges

    ionization_thresholds = spec.ionization_thresholds
    n_thresholds = length(ionization_thresholds)
    primary_transfer_matrix = zeros(n_E, n_E, n_thresholds)
    secondary_transfer_matrix = zeros(n_E, n_E, n_thresholds)

    verbose && print("Calculating energy-degradation transfer matrices for e⁻ - $(spec.name) ionizing collisions...")

    # Pre-allocate hcubature work buffers. One per thread (heap located).
    primary_bufs = [hcubature_buffer(DegradedCascadingIntegrand(0.0, 1.0, 0.0, spec.secondary_law),
                                     (0.0, 0.0), (1.0, 1.0)) for _ in 1:Threads.maxthreadid()]
    secondary_bufs = [hcubature_buffer(SecondaryCascadingIntegrand(0.0, 1.0, 0.0, spec.secondary_law),
                                       (0.0, 0.0), (1.0, 1.0)) for _ in 1:Threads.maxthreadid()]

    # Loop over ionization thresholds
    for i_threshold in n_thresholds:-1:1
        threshold = ionization_thresholds[i_threshold]

        # Find the first primary bin whose left edge is above the threshold and can
        # therefore be contribute to ionization collisions.
        i_min_primary = searchsortedfirst(E_left, threshold)

        # Loop over primary electron energy bins
        Threads.@threads :static for i_primary in i_min_primary:n_E
            E_primary_bin_min = E_edges[i_primary]      # left edge
            E_primary_bin_max = E_edges[i_primary + 1]  # right edge

            # Define the integrands for this primary bin and ionization threshold.
            degraded_integrand = DegradedCascadingIntegrand(E_primary_bin_min,
                                                           E_primary_bin_max,
                                                           threshold,
                                                           spec.secondary_law)
            secondary_integrand = SecondaryCascadingIntegrand(E_primary_bin_min,
                                                              E_primary_bin_max,
                                                              threshold,
                                                              spec.secondary_law)


            # For a fixed primary energy, the secondary/degraded boundary sits at half of
            # the excess energy. Using the lower edge of the primary bin gives the lowest
            # such boundary for this primary bin.
            E_secondary_boundary_lower = (E_primary_bin_min - threshold) / 2
            # First bin that can receive a degraded primary electron. Its left edge can be
            # ≤ E_secondary_boundary_lower as long as its right edge extends into
            # the degraded-primary range. As such, we select based on the first right edge
            # above the lowest secondary/degraded boundary for this primary bin.
            i_min_degraded = searchsortedlast(E_left, E_secondary_boundary_lower)
            # If the lowest secondary/degraded boundary lies below the grid minimum, skip.
            i_min_degraded == 0 && continue
            # Loop over degraded primary electron energy bins
            for i_degraded in i_min_degraded:(i_primary - 1)
                E_degraded_bin_min = E_edges[i_degraded]
                E_degraded_bin_max = E_edges[i_degraded + 1]
                E_degraded_lower = max(E_degraded_bin_min, E_secondary_boundary_lower)
                E_degraded_upper = min(E_degraded_bin_max, E_primary_bin_max - threshold)
                # Integrate only if limits are physical
                if E_degraded_upper > E_degraded_lower
                    result, _ = hcubature(degraded_integrand,
                                         (E_degraded_lower, 0.0),
                                         (E_degraded_upper, 1.0);
                                         buffer = primary_bufs[Threads.threadid()],
                                         )
                    primary_transfer_matrix[i_primary, i_degraded, i_threshold] = result
                end
            end


            # Using the upper edge of the primary bin gives the highest secondary/degraded
            # boundary for this primary bin.
            E_secondary_boundary_upper = (E_primary_bin_max - threshold) / 2
            # Last bin that can receive a secondary electron. Its left edge must be below
            # the highest secondary/degraded boundary for this primary bin.
            i_max_secondary = searchsortedlast(E_left, E_secondary_boundary_upper)
            # Theoritically the secondary/degraded boundary could lie above the left edge of
            # the primary bin (if we use a very coarse grid). Ensure it must be under.
            i_max_secondary = min(i_max_secondary, i_primary - 1)
            # Loop over the secondary electron energy bins
            for i_secondary in 1:i_max_secondary
                E_secondary_bin_min = E_edges[i_secondary]
                E_secondary_bin_max = E_edges[i_secondary + 1]
                E_secondary_upper = min(E_secondary_bin_max, E_secondary_boundary_upper)
                # Integrate only if limits are physical
                if E_secondary_upper > E_secondary_bin_min
                    result, _ = hcubature(secondary_integrand,
                                         (E_secondary_bin_min, 0.0),
                                         (E_secondary_upper, 1.0);
                                         buffer = secondary_bufs[Threads.threadid()],
                                         )
                    secondary_transfer_matrix[i_primary, i_secondary, i_threshold] = result
                end
            end
        end
    end

    verbose && println(" done.")
    return primary_transfer_matrix, secondary_transfer_matrix, E_edges, ionization_thresholds
end


# Load the secondary electron distribution, for a given initial primary energy index
# and ionization threshold.
function secondary_spectrum(cache::SpeciesCascadingCache, i_primary::Integer,
                            E_ionization_threshold)

    i_threshold = findmin(x -> abs(x - E_ionization_threshold), cache.ionization_thresholds)[2]
    return @view(cache.secondary_transfer_matrix[i_primary, :, i_threshold])
end

function secondary_spectrum(cache::SpeciesCascadingCache, E_primary_energy,
                            E_ionization_threshold)

    i_primary = searchsortedlast(cache.E_edges, E_primary_energy) - 1
    i_primary = clamp(i_primary, 1, size(cache.secondary_transfer_matrix, 1))
    return secondary_spectrum(cache, i_primary, E_ionization_threshold)
end


# Load the degraded primary electron distribution, for a given initial primary energy index
# and ionization threshold.
function primary_spectrum(cache::SpeciesCascadingCache, i_primary::Integer,
                          E_ionization_threshold)

    i_threshold = findmin(x -> abs(x - E_ionization_threshold), cache.ionization_thresholds)[2]
    return @view(cache.primary_transfer_matrix[i_primary, :, i_threshold])
end

function primary_spectrum(cache::SpeciesCascadingCache, E_primary_energy,
                          E_ionization_threshold)

    i_primary = searchsortedlast(cache.E_edges, E_primary_energy) - 1
    i_primary = clamp(i_primary, 1, size(cache.primary_transfer_matrix, 1))
    return primary_spectrum(cache, i_primary, E_ionization_threshold)
end





# ======================================================================================== #
#                              DISK CACHING FUNCTIONS                                      #
# ======================================================================================== #


"""
    find_cascading_file(E_edges, species_dir)

Search for a pre-computed cascading spectra file with matching energy grid.

# Arguments
- `spec::CascadingSpec` contains species name, ionization thresholds and secondary distribution law
- `E_edges`: Energy grid edges to match

# Returns
- `(file_found, filepath)`: Tuple of boolean and filepath string
"""
function find_cascading_file(spec::CascadingSpec, E_edges)
    species_dir = pkgdir(AURORA, "internal_data", "e_cascading", spec.name)
    cascading_files = readdir(species_dir)

    for filename in cascading_files
        if !endswith(filename, ".mat")
            continue
        end

        if isdir(joinpath(species_dir, filename))
            continue
        end

        filepath = joinpath(species_dir, filename)
        try
            file = matopen(filepath)

            required_keys = ["Q_primary", "Q_secondary", "E_edges", "E_ionizations"]
            if all(haskey(file, key) for key in required_keys)
                E_edges_saved = read(file, "E_edges")
                if length(E_edges) <= length(E_edges_saved) && E_edges_saved[1:length(E_edges)] == E_edges
                    close(file)
                    return (true, filepath)
                end
            end

            close(file)
        catch
            continue
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
- `(primary_transfer_matrix, secondary_transfer_matrix, E_edges, ionization_thresholds)`:
    Tuple of loaded data
"""
function load_cascading_matrices(filepath)
    println("Loading cascading-matrices from file: ", basename(filepath))
    file = matopen(filepath)
    primary_transfer_matrix = read(file, "Q_primary")
    secondary_transfer_matrix = read(file, "Q_secondary")
    ionization_thresholds = read(file, "E_ionizations")
    E_edges = read(file, "E_edges")
    close(file)

    return (primary_transfer_matrix, secondary_transfer_matrix,
            E_edges, ionization_thresholds)
end


"""
    save_cascading_matrices(primary_transfer_matrix, secondary_transfer_matrix,
                            E_edges, ionization_thresholds, species_name)

Save calculated cascading matrices to a file, located in the species-specific cascading data
directory.

# Arguments
- `primary_transfer_matrix`: Transfer matrix to save
- `secondary_transfer_matrix`: Secondary-electron transfer matrix to save
- `E_edges`: Energy grid edges used for calculations
- `ionization_thresholds`: Ionization threshold energies
- `species_name`: Name of the species (for filename)
"""
function save_cascading_matrices(primary_transfer_matrix, secondary_transfer_matrix,
                                 E_edges, ionization_thresholds, species_name)
    species_dir = pkgdir(AURORA, "internal_data", "e_cascading", "$species_name")
    filename = joinpath(species_dir,
                       string("CascadingSpec", species_name, "ionization_",
                             Dates.format(now(), "yyyymmdd-HHMMSS"),
                             ".mat"))
    file = matopen(filename, "w")
    write(file, "Q_primary", primary_transfer_matrix)
    write(file, "Q_secondary", secondary_transfer_matrix)
    write(file, "E_edges", E_edges)
    write(file, "E_ionizations", ionization_thresholds)
    close(file)
end
