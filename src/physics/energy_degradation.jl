using LoopVectorization: @tturbo

#################################################################################
#                                   Update Q                                    #
#################################################################################

function update_Q!(matrices::TransportMatrices, Ie, model::AuroraModel, t,
                   B2B_inelastic_neutrals, cascading_cache, iE, cache)

    z = model.altitude_grid.h
    ne = model.ionosphere.ne
    Te = model.ionosphere.Te
    n_neutrals_data = n_neutrals(model.ionosphere)
    σ_neutrals = model.cross_sections.σ_neutrals
    collision_levels = model.cross_sections.collision_levels
    energy_grid = model.energy_grid
    E_edges = energy_grid.E_edges
    E_centers = energy_grid.E_centers
    ΔE = energy_grid.ΔE
    Ω_beam = model.scattering.Ω_beam
    μ_center = model.pitch_angle_grid.μ_center

    n_z = length(z)
    n_μ = length(μ_center)

    Q = matrices.Q  # Extract Q for convenient access

    # e-e collisions
    if iE > 1
        loss_to_thermal_electrons!(cache.thermal_e_loss, E_centers[iE], ne, Te)
        add_thermal_electron_collisions!(@view(Q[:, :, iE - 1]), @view(Ie[:, :, iE]),
                                         cache.thermal_e_loss, ΔE[iE], n_z, n_μ)
    end

    # Get the pre-allocated ionization arrays from cache
    secondary_e_flux     = cache.secondary_e_flux
    primary_e_flux       = cache.primary_e_flux
    secondary_e_spectrum = cache.secondary_e_spectrum
    primary_e_spectrum   = cache.primary_e_spectrum

    # Loop over the neutral species
    for i in 1:length(n_neutrals_data)
        n = n_neutrals_data[i]                          # Neutral density
        σ = σ_neutrals[i]                          # Array with collision cross sections
        E_levels = collision_levels[i]            # Array with collision energy levels and number of secondary e-
        B2B_inelastic = B2B_inelastic_neutrals[i]  # Array with the probabilities of scattering from beam to beam
        species_cascading = cascading_cache[i]

        add_inelastic_collisions!(Q, Ie, z, n, σ, E_levels, B2B_inelastic, energy_grid, iE, cache)

        # Zero out the ionization arrays for this species
        fill!(secondary_e_flux[i], 0)
        fill!(primary_e_flux[i], 0)
        fill!(secondary_e_spectrum[i], 0)
        fill!(primary_e_spectrum[i], 0)

        # If the energy is too low, skip the ionization calculation (use zeros)
        idx_ionization = (E_levels[:, 2] .> 0)
        if minimum(E_levels[idx_ionization, 1]) < E_edges[iE]
            compute_ionization_flux!(secondary_e_flux[i], primary_e_flux[i],
                                     n, Ie, t, z, μ_center, Ω_beam, iE, cache, i)
            compute_ionization_spectra!(secondary_e_spectrum[i], primary_e_spectrum[i],
                                        σ, E_levels, species_cascading, energy_grid, iE)
        end
    end
    # If there is no ionization to add (everything is zero), skip the update of Q
    # Mmh the compiler seems to be smart enough to skip the update of Q anyway when
    # everything is zero. But let's keep this if statement just in case.
    has_ionization = any(!iszero, secondary_e_spectrum) || any(!iszero, primary_e_spectrum)
    if has_ionization
        add_ionization_collisions!(Q, iE,
                                   secondary_e_flux, primary_e_flux,
                                   secondary_e_spectrum, primary_e_spectrum)
    end
end


function add_thermal_electron_collisions!(Q_slice, Ie_slice, thermal_e_loss, ΔE, n_z, n_μ)
    @tturbo for i_t in axes(Ie_slice, 2)
        for i_μ in 1:n_μ
            for i_z in 1:n_z
                row = (i_μ - 1) * n_z + i_z
                Q_slice[row, i_t] += thermal_e_loss[i_z] * Ie_slice[row, i_t] / ΔE
            end
        end
    end

    return nothing
end




#################################################################################
#                               Inelastic collisions                            #
#################################################################################

"""
    calculate_scattered_flux!(result, B2B_inelastic, n, Ie_slice)

Calculate the flux of electrons after pitch-angle scattering by inelastic collisions.

# Physics
The scattered flux at each altitude and angle is computed as:
    result[z, μ₁, t] = Σ_μ₂ n(z) x P(μ₁←μ₂) x Ie[z, μ₂, t]

where:
- `n(z)` is the neutral density at altitude z
- `P(μ₁←μ₂)` is the probability of scattering from pitch angle μ₂ to μ₁ (from B2B_inelastic)
- `Ie[z, μ₂, t]` is the incident electron flux before scattering

# Implementation
This exploits the block-diagonal structure of the full scattering matrix to avoid storing
and accessing a large sparse matrix. Instead, it computes the multiplication directly from
the small scattering probability matrix and density profile.

# Arguments
- `result`: Output array (n_z x n_μ, n_t) - scattered electron flux
- `B2B_inelastic`: Scattering probability matrix (n_μ, n_μ) - pitch angle redistribution
- `n`: Neutral density profile (n_z,) - altitude-dependent density [m⁻³]
- `Ie_slice`: Incident electron flux at the current energy (n_z x n_μ, n_t) - flux before scattering
"""
function calculate_scattered_flux!(result, B2B_inelastic, n, Ie_slice)
    # Zero out result first
    # fill!(result, 0.0) # actually we don't need this because all values will get updated

    n_z = length(n)
    n_μ = size(B2B_inelastic, 2)
    n_t = size(Ie_slice, 2)

    @tturbo for it in 1:n_t
        for i1 in 1:n_μ
            for iz in 1:n_z
                row = (i1 - 1) * n_z + iz
                tmp = 0.0
                for i2 in 1:n_μ
                    col = (i2 - 1) * n_z + iz
                    tmp += n[iz] * B2B_inelastic[i1, i2] * Ie_slice[col, it]
                end
                result[row, it] = tmp
            end
        end
    end

    return nothing
end

function add_inelastic_collisions!(Q, Ie, z, n, σ, E_levels, B2B_inelastic, energy_grid::EnergyGrid, iE, cache)
    E_edges = energy_grid.E_edges
    ΔE = energy_grid.ΔE
    # Use cached matrices to avoid allocations
    Ie_scatter = cache.Ie_scatter

    # Calculate the flux of electrons after pitch-angle scattering by inelastic collisions
    # This is computed ONCE for all energy levels of this species at this energy
    calculate_scattered_flux!(Ie_scatter, B2B_inelastic, n, @view(Ie[:, :, iE]))

    # Loop over the energy levels of the collisions with the i-th neutral species
    for i_level in axes(E_levels, 1)[2:end]
        if E_levels[i_level, 2] <= 0  # these collisions should not produce secondary e-
            # Calculate the degradation factor combining:
            # 1) Cross-section σ[i_level, iE] - collision probability
            # 2) Energy loss correction min(1, E_loss/ΔE[iE]) - accounts for when the energy
            #    loss is smaller than the bin width (prevents over-depleting the current bin)
            factor = σ[i_level, iE] * min(1, E_levels[i_level, 1] / ΔE[iE])

            # Find the energy bins where electrons will end up after losing E_levels[i_level, 1] eV
            E_loss = E_levels[i_level, 1]
            E_min = E_edges[iE] - E_loss           # Minimum energy after collision
            E_max = E_edges[iE+1] - E_loss         # Maximum energy after collision
            # Find indices where degraded electrons end up
            i_min = searchsortedfirst(@view(E_edges[2:end]), E_min)  # First bin with upper edge > E_min
            i_max = searchsortedlast(@view(E_edges[1:end-1]), E_max) # Last bin with lower edge < E_max
            i_degrade = i_min:i_max
            partition_fraction = zeros(length(i_degrade)) # initialise

            if !isempty(i_degrade) && i_degrade[1] < iE
                # Distribute the degrading e- between those bins
                partition_fraction[1] = min(1, (E_edges[i_degrade[1]+1] -
                                                E_edges[iE] + E_levels[i_level, 1]) / ΔE[iE])
                if length(i_degrade) > 2
                    partition_fraction[2:end-1] = min.(1, ΔE[i_degrade[2:end-1]] / ΔE[iE])
                end
                partition_fraction[end] = min(1, (E_edges[iE+1] - E_edges[i_degrade[end]] -
                                                    E_levels[i_level, 1]) / ΔE[iE])
                if i_degrade[end] == iE
                    partition_fraction[end] = 0
                end

                # Normalize partition fractions to sum to 1
                partition_fraction = partition_fraction / sum(partition_fraction)

                # Add the degraded electron flux to Q
                # Q[z,t,E'] += Ie_scatter[z,t] x partition_fraction[E'] x σ x min(1, E_loss/dE)
                @tturbo for i_u in eachindex(findall(x -> x != 0, partition_fraction))
                    for j in axes(Q, 2)
                        for k in axes(Q, 1)
                            Q[k, j, i_degrade[i_u]] += Ie_scatter[k, j] * partition_fraction[i_u] * factor
                        end
                    end
                end
            end
        end
    end
end




#################################################################################
#                               Ionization collisions                           #
#################################################################################

#=
The three functions below compute the ionization contribution to Q more efficiently
than a straightforward per-channel loop would.

The key insight is that for a given neutral species, the only quantities that vary
across ionization channels are the cross-section σ[channel, iE] and the energy
spectra for secondaries and degraded primaries. The spatial/angular structure of the
electron flux is the same for all channels of the same species, so we factorize the
calculation into two independent parts.

PART 1 — Spatial/angular flux (computed once per species with `compute_ionization_flux!`)
  - `secondary_e_flux[i_s]` (shape: n_z·n_μ × n_t): incident electron flux redistributed
    isotropically over all pitch angles, weighted by the neutral density.
        secondary_e_flux = n × Ie[:, :, iE]    (uniformly spread over μ via solid angles)
  - `primary_e_flux[i_s]` (shape: n_z·n_μ × n_t): incident electron flux staying in
    its original pitch-angle beam, weighted by the neutral density.
        primary_e_flux = n × Ie[:, :, iE]    (same pitch angle as the incident electron)

PART 2 — Energy spectrum weights (computed once per species with `compute_ionization_spectra!`)
  - `secondary_e_spectrum[i_s]` (shape: n_E): cross-section-weighted sum of the secondary
    electron energy spectra over all ionization channels of the species:
        secondary_e_spectrum = ∑_channels  secondary_e_spectra[channel] × σ[channel, iE]
  - `primary_e_spectrum[i_s]` (shape: n_E): same but for the degraded primary electrons:
        primary_e_spectrum = ∑_channels  primary_e_spectra[channel] × σ[channel, iE]

The factorization in Part 2 is the key optimization: all ionization channels of a species
are collapsed into a single spectrum vector, so the expensive (n_E × n_z × n_μ × n_t)
accumulation into Q only needs to run once per species instead of once per channel.

The final accumulation into Q is done in `add_ionization_collisions!`:
```julia
for i_s in eachindex(secondary_e_flux)
    for iI in 1:(iE - 1)
        for j in axes(Q, 2)
            for k in axes(Q, 1)
                Q[k, j, iI] += secondary_e_flux[i_s][k, j] * secondary_e_spectrum[i_s][iI] +  # secondaries of species i_s
                               primary_e_flux[i_s][k, j]   * primary_e_spectrum[i_s][iI]      # degraded primaries of species i_s
            end
        end
    end
end
```
=#
function compute_ionization_flux!(secondary_e_flux, primary_e_flux,
                                  n, Ie, t, z, μ_center, Ω_beam, iE,
                                  cache, i_species)
    source_sum = cache.ionization_source_sum

    n_z = length(z)
    n_μ = length(μ_center)
    n_t = size(Ie, 2)
    inv_sum_Ω_beam = inv(sum(Ω_beam))

    # PRIMARY ELECTRONS: simply the incident flux in each pitch-angle beam, weighted by the
    # neutral density. The primary electrons continue propagating in the same direction (beam)
    # after the collision.
    @tturbo for it in 1:n_t
        for i_μ in 1:n_μ
            for iz in 1:n_z
                row = (i_μ - 1) * n_z + iz
                primary_e_flux[row, it] = n[iz] * Ie[row, it, iE]
            end
        end
    end

    # SECONDARY ELECTRONS: we first compute the total incident flux (summed over pitch-angle beams)
    @tturbo for it in 1:n_t
        for iz in 1:n_z
            source_total = 0.0
            for i_μ in 1:n_μ
                row = (i_μ - 1) * n_z + iz
                source_total += Ie[row, it, iE]
            end
            source_sum[iz, it] = source_total
        end
    end

    # SECONDARY ELECTRONS: which we then redistribute isotropically over all pitch angles,
    # using the solid-angle weight fractions Ω_beam[μᵢ] / sum(Ω_beam). We also weight by the
    # neutral density.
    @tturbo for it in 1:n_t
        for i_μ in 1:n_μ
            beam_weight = Ω_beam[i_μ] * inv_sum_Ω_beam
            for iz in 1:n_z
                row = (i_μ - 1) * n_z + iz
                secondary_e_flux[row, it] = n[iz] * source_sum[iz, it] * beam_weight
            end
        end
    end

    return nothing
end

function compute_ionization_spectra!(secondary_e_spectrum, primary_e_spectrum,
                                     σ, E_levels, species_cascading,
                                     energy_grid::EnergyGrid, iE)
    # Loop through the ionization channels (E_levels[:,2] > 0) for the current neutral species.
    # For each channel, compute the energy spectra for secondaries and degraded primaries, then
    # accumulate them weighted by the channel's cross-section. This factorizes all channels into
    # a single spectrum vector per species.
    for i_level in axes(E_levels, 1)[2:end]
        if E_levels[i_level, 2] > 0    # ionizing collision → produces secondary electrons
            E_loss = E_levels[i_level, 1]
            # Retrieve precomputed, bin-integrated spectra from the cascading cache.
            secondary_e_spectra = secondary_spectrum(species_cascading, iE, E_loss)
            primary_e_spectra = primary_spectrum(species_cascading, iE, E_loss)

            if sum(secondary_e_spectra) > 0
                secondary_e_spectra = secondary_e_spectra ./ sum(secondary_e_spectra) # normalize sum to 1
                secondary_e_spectra = E_levels[i_level, 2] .* secondary_e_spectra  # scale by number of secondaries
            end
            if sum(primary_e_spectra) > 0
                primary_e_spectra = primary_e_spectra ./ sum(primary_e_spectra) # normalize sum to 1
            end

            # Accumulate into the species-level spectrum, weighted by cross-section
            secondary_e_spectrum .+= secondary_e_spectra .* σ[i_level, iE]
            primary_e_spectrum   .+= primary_e_spectra   .* σ[i_level, iE]
        end
    end
end

@generated function add_ionization_collisions!(Q, iE,
                                               secondary_e_flux::NTuple{N, <:AbstractMatrix},
                                               primary_e_flux::NTuple{N, <:AbstractMatrix},
                                               secondary_e_spectrum::NTuple{N, <:AbstractVector},
                                               primary_e_spectrum::NTuple{N, <:AbstractVector}) where {N}

    # Unroll the loop over species to accumulate in Q only once and reduce memory access.
    terms = [:(secondary_e_flux[$i_s][k, j] * secondary_e_spectrum[$i_s][iI] +
               primary_e_flux[$i_s][k, j] * primary_e_spectrum[$i_s][iI])
             for i_s in 1:N]
    rhs = reduce((a, b) -> :($a + $b), terms)

    quote
        @tturbo for iI in 1:(iE - 1)
            for j in axes(Q, 2)
                for k in axes(Q, 1)
                    Q[k, j, iI] += $rhs
                end
            end
        end
    end
end
