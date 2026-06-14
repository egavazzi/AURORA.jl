using LoopVectorization: @tturbo

#################################################################################
#                                   Update Q                                    #
#################################################################################

function update_Q!(matrices::TransportMatrices, Ie, model::AuroraModel, t,
                   B2B_inelastic_neutrals, iE, cache)

    z = model.altitude_grid.h
    ne = model.ionosphere.ne
    Te = model.ionosphere.Te
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

    # Get the pre-allocated degradation arrays from cache
    Ie_scatter           = cache.Ie_scatter
    inelastic_weight     = cache.inelastic_weight
    secondary_e_flux     = cache.secondary_e_flux
    primary_e_flux       = cache.primary_e_flux
    secondary_e_spectrum = cache.secondary_e_spectrum
    primary_e_spectrum   = cache.primary_e_spectrum

    # Loop over the neutral species
    for (i, sp) in enumerate(model.species)
        n = sp.density
        σ = sp.cross_sections
        E_levels = sp.excitation_levels
        B2B_inelastic = B2B_inelastic_neutrals[i]
        species_cascading = sp.cascading_data

        # ── Inelastic (non-ionizing) collisions ──
        # Compute the scattered flux and the per-target-bin degradation weights, but
        # do NOT accumulate into Q here. The accumulation is fused below across all
        # species and channels so that Q is swept only once.
        calculate_scattered_flux!(Ie_scatter[i], B2B_inelastic, n, @view(Ie[:, :, iE]))
        compute_inelastic_weights!(inelastic_weight[i], σ, E_levels, energy_grid, iE)

        # Zero out the ionization arrays for this species
        fill!(secondary_e_flux[i], 0)
        fill!(primary_e_flux[i], 0)
        fill!(secondary_e_spectrum[i], 0)
        fill!(primary_e_spectrum[i], 0)

        # If the energy is too low, skip the ionization calculation (use zeros).
        # Find the minimum ionization-channel threshold without allocating a mask/index array.
        min_ionization_E = Inf
        for i_level in axes(E_levels, 1)
            if E_levels[i_level, 2] > 0
                min_ionization_E = min(min_ionization_E, E_levels[i_level, 1])
            end
        end
        if min_ionization_E < E_edges[iE]
            compute_ionization_flux!(secondary_e_flux[i], primary_e_flux[i],
                                     n, Ie, z, μ_center, Ω_beam, iE, cache)
            compute_ionization_spectra!(secondary_e_spectrum[i], primary_e_spectrum[i],
                                        σ, E_levels, species_cascading, iE)
        end
    end
    # Fused accumulation into Q for ALL degradation channels (inelastic losses +
    # ionization secondaries + degraded primaries) of ALL species. This touches the
    # lower-energy block of Q exactly once per energy step instead of once per
    # species/excitation-level/target-bin. Skipped entirely when there is nothing to
    # degrade (lowest energy bin, or all contributions zero).
    has_degradation = any(!iszero, inelastic_weight) ||
                      any(!iszero, secondary_e_spectrum) ||
                      any(!iszero, primary_e_spectrum)
    if iE > 1 && has_degradation
        add_degradation_collisions!(Q, iE,
                                    Ie_scatter, inelastic_weight,
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

# Partition fraction for position `u` (1-based) within the target bin range `i_degrade`,
# replicating the original first/middle/last bin handling exactly. The last position is
# evaluated first so that a single-bin range (u == 1 == n_deg) uses the "last" formula,
# matching the original code where partition_fraction[end] overwrote partition_fraction[1].
@inline function inelastic_partition_fraction(u, n_deg, i_degrade, iE, E_loss, E_edges, ΔE)
    if u == n_deg
        idx_end = i_degrade[n_deg]
        idx_end == iE && return 0.0
        return min(1.0, (E_edges[iE + 1] - E_edges[idx_end] - E_loss) / ΔE[iE])
    elseif u == 1
        return min(1.0, (E_edges[i_degrade[1] + 1] - E_edges[iE] + E_loss) / ΔE[iE])
    else
        return min(1.0, ΔE[i_degrade[u]] / ΔE[iE])
    end
end

"""
    compute_inelastic_weights!(weight, σ, E_levels, energy_grid, iE)

Compute, for one neutral species, the per-target-bin energy-degradation weights from
non-ionizing inelastic collisions of an electron at energy index `iE`.

`weight[iE_degrade]` accumulates, over every non-ionizing excitation channel,
```
    partition_fraction × σ[level, iE] × min(1, E_loss / ΔE[iE])
```
where `partition_fraction` distributes the degraded electrons over the lower-energy bins
they fall into. The weights are consumed once by the fused `add_degradation_collisions!`
together with the species' scattered flux `Ie_scatter`, so that
```
    Q[:, :, iE_degrade] += Ie_scatter[:, :] × weight[iE_degrade]
```
This replaces the original per-channel accumulation directly into Q.

# Arguments
- `weight`: output vector (n_E,), overwritten with the degradation weights
- `σ`: cross-sections for this species (n_levels × n_E)
- `E_levels`: excitation-level table (column 1 = energy loss, column 2 = #secondaries)
- `energy_grid::EnergyGrid`
- `iE`: current energy index
"""
function compute_inelastic_weights!(weight, σ, E_levels, energy_grid::EnergyGrid, iE)
    E_edges = energy_grid.E_edges
    ΔE = energy_grid.ΔE

    fill!(weight, 0.0)

    # Loop over the energy levels of the collisions with this neutral species
    for i_level in axes(E_levels, 1)[2:end]
        E_levels[i_level, 2] <= 0 || continue  # only collisions that do not produce secondary e-

        # Degradation factor combining:
        # 1) Cross-section σ[i_level, iE] - collision probability
        # 2) Energy loss correction min(1, E_loss/ΔE[iE]) - accounts for when the energy
        #    loss is smaller than the bin width (prevents over-depleting the current bin)
        E_loss = E_levels[i_level, 1]
        factor = σ[i_level, iE] * min(1.0, E_loss / ΔE[iE])

        # Find the energy bins where electrons will end up after losing E_loss eV
        E_min = E_edges[iE] - E_loss           # Minimum energy after collision
        E_max = E_edges[iE + 1] - E_loss       # Maximum energy after collision
        i_min = searchsortedfirst(@view(E_edges[2:end]), E_min)   # First bin with upper edge > E_min
        i_max = searchsortedlast(@view(E_edges[1:end - 1]), E_max) # Last bin with lower edge < E_max
        i_degrade = i_min:i_max
        (isempty(i_degrade) || i_degrade[1] >= iE) && continue

        n_deg = length(i_degrade)

        # Normalize the partition fractions to sum to 1, then accumulate into `weight`.
        sum_pf = 0.0
        for u in 1:n_deg
            sum_pf += inelastic_partition_fraction(u, n_deg, i_degrade, iE, E_loss, E_edges, ΔE)
        end
        sum_pf > 0 || continue
        norm_factor = factor / sum_pf
        for u in 1:n_deg
            pf = inelastic_partition_fraction(u, n_deg, i_degrade, iE, E_loss, E_edges, ΔE)
            weight[i_degrade[u]] += pf * norm_factor
        end
    end

    return nothing
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

The final accumulation into Q is done in `add_degradation_collisions!`, which fuses these
ionization terms together with the inelastic degradation term (`Ie_scatter × inelastic_weight`,
see `compute_inelastic_weights!`) into a single pass over Q:
```julia
for i_s in eachindex(secondary_e_flux)
    for iI in 1:(iE - 1)
        for j in axes(Q, 2)
            for k in axes(Q, 1)
                Q[k, j, iI] += Ie_scatter[i_s][k, j]     * inelastic_weight[i_s][iI] +     # inelastic degradation of species i_s
                               secondary_e_flux[i_s][k, j] * secondary_e_spectrum[i_s][iI] + # secondaries of species i_s
                               primary_e_flux[i_s][k, j]   * primary_e_spectrum[i_s][iI]      # degraded primaries of species i_s
            end
        end
    end
end
```
=#
function compute_ionization_flux!(secondary_e_flux, primary_e_flux,
                                  n, Ie, z, μ_center, Ω_beam, iE,
                                  cache)
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
                                     σ, E_levels, species_cascading, iE)
    # Loop through the ionization channels (E_levels[:,2] > 0) for the current neutral species.
    # For each channel, compute the energy spectra for secondaries and degraded primaries, then
    # accumulate them weighted by the channel's cross-section. This factorizes all channels into
    # a single spectrum vector per species.
    for i_level in axes(E_levels, 1)[2:end]
        if E_levels[i_level, 2] > 0    # ionizing collision → produces secondary electrons
            E_loss = E_levels[i_level, 1]
            n_secondary = E_levels[i_level, 2]
            σ_level = σ[i_level, iE]
            # Retrieve precomputed, bin-integrated spectra from the cascading cache.
            secondary_e_spectra = secondary_spectrum(species_cascading, iE, E_loss)
            primary_e_spectra = primary_spectrum(species_cascading, iE, E_loss)

            sum_secondary = sum(secondary_e_spectra)    # for normalization
            sum_primary = sum(primary_e_spectra)        # for normalization
            if sum_secondary > 0
                # scale by cross-section, spectra normalization, and number of secondaries
                secondary_scale = σ_level * n_secondary / sum_secondary
                secondary_e_spectrum .+= secondary_e_spectra .* secondary_scale
            end
            if sum_primary > 0
                # scale by cross-section and spectra normalization
                primary_scale = σ_level / sum_primary
                primary_e_spectrum .+= primary_e_spectra .* primary_scale
            end
        end
    end
end

@generated function add_degradation_collisions!(Q, iE,
                                                Ie_scatter::NTuple{N, <:AbstractMatrix},
                                                inelastic_weight::NTuple{N, <:AbstractVector},
                                                secondary_e_flux::NTuple{N, <:AbstractMatrix},
                                                primary_e_flux::NTuple{N, <:AbstractMatrix},
                                                secondary_e_spectrum::NTuple{N, <:AbstractVector},
                                                primary_e_spectrum::NTuple{N, <:AbstractVector}) where {N}

    # Unroll the loop over species so that Q is accumulated in a single pass, fusing all
    # three degradation channels (inelastic losses + ionization secondaries + degraded
    # primaries) and reducing memory traffic on the large Q array.
    terms = [:(Ie_scatter[$i_s][k, j] * inelastic_weight[$i_s][iI] +
               secondary_e_flux[$i_s][k, j] * secondary_e_spectrum[$i_s][iI] +
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
