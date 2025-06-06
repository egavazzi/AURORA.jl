using LoopVectorization: @tturbo
using SparseArrays: sparse!



#################################################################################
#                                   Update Q                                    #
#################################################################################

function demo_update_Q!(Q, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals,
                   B2B_inelastic_neutrals, cascading_neutrals, E, dE, iE, BeamWeight,
                   μ_center, cache)

    # e-e collisions
    @views if iE > 1
        Q[:, :, iE - 1] .+= repeat(loss_to_thermal_electrons(E[iE], ne, Te) / dE[iE],
                                outer = (length(μ_center), length(t))) .* Ie[:, :, iE];
    end

    # Loop over the species
    for i in 1:length(n_neutrals)
        n = n_neutrals[i];                          # Neutral density
        σ = σ_neutrals[i];                          # Array with collision cross sections
        E_levels = E_levels_neutrals[i];            # Array with collision enery levels and number of secondary e-
        B2B_inelastic = B2B_inelastic_neutrals[i];  # Array with the probablities of scattering from beam to beam
        cascading = cascading_neutrals[i];          # Cascading function for the current i-th species

        add_inelastic_collisions!(Q, Ie, h_atm, n, σ, E_levels, B2B_inelastic, E, dE, iE, cache)
        add_ionization_collisions!(Q, Ie, h_atm, t, n, σ, E_levels, cascading, E, dE, iE, BeamWeight, μ_center)
    end
end



function update_Q!(Q, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals,
                      B2B_inelastic_neutrals, cascading_neutrals, E, dE, iE, BeamWeight,
                      μ_center, cache)

    # e-e collisions
    @views if iE > 1
        Q[:, :, iE - 1] .+= repeat(loss_to_thermal_electrons(E[iE], ne, Te) / dE[iE],
                                   outer = (length(μ_center), length(t))) .* Ie[:, :, iE]
    end

    # TODO for later: move this outside of the energy loop and then just fill with zeros
    Ionization_fragment_1 = [zeros(size(Ie, 1), size(Ie, 2)) for _ in 1:length(n_neutrals)]
    Ionizing_fragment_1 = [zeros(size(Ie, 1), size(Ie, 2)) for _ in 1:length(n_neutrals)]
    Ionization_fragment_2 = [zeros(length(E)) for _ in 1:length(n_neutrals)]
    Ionizing_fragment_2 = [zeros(length(E)) for _ in 1:length(n_neutrals)]

    # Loop over the neutral species
    for i in 1:length(n_neutrals)
        n = n_neutrals[i]                          # Neutral density
        σ = σ_neutrals[i]                          # Array with collision cross sections
        E_levels = E_levels_neutrals[i]            # Array with collision enery levels and number of secondary e-
        B2B_inelastic = B2B_inelastic_neutrals[i]  # Array with the probablities of scattering from beam to beam
        cascading = cascading_neutrals[i]          # Cascading function for the current i-th species

        add_inelastic_collisions!(Q, Ie, h_atm, n, σ, E_levels, B2B_inelastic, E, dE, iE, cache)

        # TODO: that's if I will move things outside of the energy loop
        # fill!(Ionization_fragment_1[i], 0)
        # fill!(Ionizing_fragment_1[i], 0)
        # fill!(Ionization_fragment_2[i], 0)
        # fill!(Ionizing_fragment_2[i], 0)

        # If the energy is too low, skip the ionization calculation (use zeros)
        idx_ionization = (E_levels[:, 2] .> 0)
        if minimum(E_levels[idx_ionization, 1]) < E[iE]
            prepare_first_ionization_fragment!(Ionization_fragment_1[i], Ionizing_fragment_1[i],
                                    n, Ie, t, h_atm, μ_center, BeamWeight, iE)
            prepare_second_ionization_fragment!(Ionization_fragment_2[i], Ionizing_fragment_2[i],
                                    σ, E_levels, cascading, E, dE, iE)
        end
    end
    # If there is no ionization to add (everything is zero), skip the update of Q
    # Mmh the compiler seems to be smart enough to skip the update of Q anyway when
    # everything is zero. But let's keep this if statement just in case.
    if !(iszero(Ionization_fragment_2) && iszero(Ionizing_fragment_2))
        add_ionization_fragments!(Q, iE,
                                  Ionization_fragment_1, Ionizing_fragment_1,
                                  Ionization_fragment_2, Ionizing_fragment_2)
    end
end





#################################################################################
#                               Inelastic collisions                            #
#################################################################################

# legacy demo function
function make_big_B2B_matrix(B2B_inelastic, n, h_atm)
    idx1 = Vector{Int}(undef, 0)
    idx2 = Vector{Int}(undef, 0)
    aB2B = Vector{Float64}(undef, 0)
    for i1 in axes(B2B_inelastic, 1)
        for i2 in axes(B2B_inelastic, 2)
            append!(idx1, (i1 - 1) * length(h_atm) .+ (1:length(h_atm)))
            append!(idx2, (i2 - 1) * length(h_atm) .+ (1:length(h_atm)))
            append!(aB2B, n * B2B_inelastic[i1, i2])
        end
    end
    AB2B = sparse!(idx1, idx2, aB2B)
    return AB2B
end

# New function, this is run one time (typically at first iteration of the energy loop) to
# fully build the sparse B2B matrix
function build_big_B2B(B2B_inelastic, n, h_atm)
    n_elements = size(B2B_inelastic, 1) * size(B2B_inelastic, 2) * length(h_atm)
    idx1 = Vector{Int}(undef, n_elements)
    idx2 = Vector{Int}(undef, n_elements)
    aB2B = Vector{Float64}(undef, n_elements)
    counter = 1
    for i1 in axes(B2B_inelastic, 1)
        for i2 in axes(B2B_inelastic, 2)
            idx_range = counter:(counter + length(h_atm) - 1)
            idx1[idx_range] .= (i1 - 1) * length(h_atm) .+ (1:length(h_atm))
            idx2[idx_range] .= (i2 - 1) * length(h_atm) .+ (1:length(h_atm))
            aB2B[idx_range] .= n * B2B_inelastic[i1, i2]
            counter += length(h_atm)
        end
    end

    return sparse!(idx1, idx2, aB2B)
end

# New function, this is run every time we need to update the values in the B2B matrix (all
# other iterations of the energy loop)
# This was made possible after some trial and error to find out in what order the elements
# of the B2B sparse matrix are stored in memory. Here we update these values in place.
# What would'nt we do for better performance eh?
function update_big_B2B!(AB2B, B2B_inelastic, n)
    counter = 1
    n_μ = size(B2B_inelastic, 2)
    @views for i1 in axes(B2B_inelastic, 2)
        for i2 in eachindex(n)
            idx_range = counter:(counter + n_μ - 1)
            AB2B.nzval[idx_range] .= n[i2] .* B2B_inelastic[:, i1]
            counter += n_μ
        end
    end
    return nothing
end



function add_inelastic_collisions!(Q, Ie, h_atm, n, σ, E_levels, B2B_inelastic, E, dE, iE, cache)
    # Initialize
    Ie_degraded = Matrix{Float64}(undef, size(Ie, 1), size(Ie, 2))

    # Multiply each element of B2B with n (density vector) and resize to get a matrix that
    # can be multiplied with Ie
    if iE == length(E)
        # build the matrix the first time
        cache.AB2B = build_big_B2B(B2B_inelastic, n, h_atm)
    else
        # then reuse the sparse structure and just update the values
        update_big_B2B!(cache.AB2B, B2B_inelastic, n)
    end

    # Then calculate the flux of electrons that scatter
    Ie_scatter = cache.AB2B * @view(Ie[:, :, iE])

    # Loop over the energy levels of the collisions with the i-th neutral species
    for i_level in axes(E_levels, 1)[2:end]
        if E_levels[i_level, 2] <= 0  # these collisions should not produce secondary e-
            # The flux of e- degraded from energy bin [E[iE], E[iE] + dE[iE]] to any lower
            # energy bin by excitation of the E_levels[i_level] state of the current
            # species. The second factor corrects for the case where the energy loss is
            # smaller than the width in energy of the energy bin. That is, when dE[iE] >
            # E_levels[i_level,1], only the fraction E_levels[i_level,1] / dE[iE] is lost
            # from the energy bin [E[iE], E[iE] + dE[iE]].

            Ie_degraded .= (σ[i_level, iE] * min(1, E_levels[i_level, 1] / dE[iE])) .* Ie_scatter

            # Find the energy bins where the e- in the current energy bin will degrade when
            # losing E_levels[i_level, 1] eV
            i_degrade = intersect(findall(x -> x > E[iE] - E_levels[i_level, 1], E + dE),      # find lowest bin
                                  findall(x -> x < E[iE] + dE[iE] - E_levels[i_level, 1], E))  # find highest bin

            partition_fraction = zeros(size(i_degrade)) # initialise
            if !isempty(i_degrade) && i_degrade[1] < iE
                # Distribute the degrading e- between those bins
                partition_fraction[1] = min(1, (E[i_degrade[1]] .+ dE[i_degrade[1]] .-
                                                E[iE] .+ E_levels[i_level, 1]) / dE[iE])
                if length(i_degrade) > 2
                    partition_fraction[2:end-1] = min.(1, dE[i_degrade[2:end-1]] / dE[iE])
                end
                partition_fraction[end] = min(1, (E[iE] .+ dE[iE] .- E[i_degrade[end]] .-
                                                    E_levels[i_level, 1]) / dE[iE])
                if i_degrade[end] == iE
                    partition_fraction[end] = 0
                end

                # normalise
                partition_fraction = partition_fraction / sum(partition_fraction)

                # and finally calculate the flux of degrading e-
                @tturbo for i_u in eachindex(findall(x -> x != 0, partition_fraction))
                    for j in axes(Q, 2)
                        for k in axes(Q, 1)
                            Q[k, j, i_degrade[i_u]] +=  max(0, Ie_degraded[k, j]) * partition_fraction[i_u]
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
This is kind of the demo function for the ionization (legacy).
For a given energy E and a neutral species n, this function does the following:

1. Loop over the ionization collisions of the current neutral species n
2. For each ionization collision, calculate the ionization matrix and ionizing matrix
3. For each ionization collision, calculate the secondary_e_spectra vector and primary_e_spectra vector
4. Add everything to Q

Let's review the different elements:

## Ionization matrix, shape (n_μ x n_z, n_μ x n_z)
It is the flux of secondary electrons. It is calculated as the product of the density of the
neutral species, the cross-section of the current ionization collision, and the electron
flux at the current energy. It is also redristibuted isotropically over the different angles.
A simplified calculation looks like this:
    Ionization_matrix = n * σ * Ie[:, :, iE]    (no isotropic redistribution here)

## Ionizing matrix, shape (n_μ x n_z, n_μ x n_z)
It is the flux of primary electrons. It is calculated as the product of the density of the
neutral species, the cross-section of the current ionization collision, and the electron
flux at the current energy. The difference with the flux of secondary electron (Ionization
matrix) is that it is not redistributed (the primary electrons continue to propagate with
the same pitch-angle as before).
A simplified calculation looks like this:
    Ionizing_matrix = n * σ * Ie[:, :, iE]

## secondary_e_spectra and primary_e_spectra vectors, both of shape (nE)
These are calculated using the very fancy `cascading` function. The function calculates how
the secondary and degraded primary electrons will be distributed in energy after the collision.
They are scaled so that sum(secondary_e_spectra) = number_of_secondaries (1 or 2 depending
on the collision) and sum(primary_e_spectra) = 1.


The (Ionization matrix + secondary_e_spectra) and the (Ionizing matrix + primary_e_spectra)
form two pairs. The first is used to calculate the final secondary electrons that will be
created at lower energies, and the second one to calculate the final primary electrons that
will be degraded in energy.

These pairs are added to Q using a for loop of the following form

```julia
for iI in 1:(iE - 1)
    for j in axes(Q, 2)
        for k in axes(Q, 1)
            Q[k, j, iI] += Ionization[k, j] * secondary_e_spectra[iI] +
                            Ionizing[k, j] * primary_e_spectra[iI]
        end
    end
end
```
=#
function add_ionization_collisions!(Q, Ie, h_atm, t, n, σ, E_levels, cascading, E, dE, iE,
                                    BeamWeight, μ_center)

    Ionization = zeros(size(Ie, 1), size(Ie, 2))
    Ionizing = Matrix{Float64}(undef, size(Ie, 1), size(Ie, 2))
    n_repeated_over_μt = repeat(n, length(μ_center), length(t))
    n_repeated_over_t = repeat(n, 1, length(t))

    for i_level in axes(E_levels, 1)[2:end]
        if E_levels[i_level, 2] > 0    # these collisions should produce secondary e-
            # Find the energy bins where the e- in the current energy bin will degrade when
            # losing E_levels[i_level, 1] eV
            i_degrade = intersect(findall(x -> x > E[iE] - E_levels[i_level, 1], E + dE),     # find lowest bin
                                  findall(x -> x < E[iE] + dE[iE] - E_levels[i_level,1], E))  # find highest bin

            if !isempty(i_degrade) && i_degrade[1] < iE
                # ISOTROPIC SECONDARY ELECTRONS

                Ionizing .= n_repeated_over_μt .* (σ[i_level, iE] .* @view(Ie[:, :, iE]));
                fill!(Ionization, 0)
                @views for i_μ1 in eachindex(μ_center)
                    for i_μ2 in eachindex(μ_center)
                        Ionization[(i_μ1 - 1) * length(h_atm) .+ (1:length(h_atm)), :] .+=
                            max.(0, n_repeated_over_t .*
                                (σ[i_level, iE] .*
                                Ie[(i_μ2 - 1) * length(h_atm) .+ (1:length(h_atm)), :, iE]) .*
                                BeamWeight[i_μ1] ./ sum(BeamWeight))
                    end
                end

                # Calculate the spectra of the secondary e-
                secondary_e_spectra = cascading(E, E[iE], E_levels[i_level, 1], "s");
                # We use the average energy of the e- in the current energy bin
                secondary_e_spectra = (secondary_e_spectra .+ secondary_e_spectra[[2:end; end]]) .* dE / 2

                # Calculate the distribution of the ionizing (= primary) e-, that have lost the
                # corresponding amount of energy
                primary_e_spectra = cascading(E, E[iE], E_levels[i_level, 1], "c");

                if sum(secondary_e_spectra) > 0
                    # normalise
                    secondary_e_spectra = secondary_e_spectra ./ sum(secondary_e_spectra)
                    # and multiply with the number of e- created
                    secondary_e_spectra = E_levels[i_level, 2] .* secondary_e_spectra
                end
                if sum(primary_e_spectra) > 0
                    # normalise
                    primary_e_spectra = primary_e_spectra ./ sum(primary_e_spectra)
                end

                # and finally add this to the flux of degrading e-
                @tturbo for iI in 1:(iE - 1)
                    for j in axes(Q, 2)
                        for k in axes(Q, 1)
                            Q[k, j, iI] += Ionization[k, j] * secondary_e_spectra[iI] +
                                           Ionizing[k, j] * primary_e_spectra[iI]
                        end
                    end
                end
            end
        end
    end
end



#=
The three following functions do the same thing as the function above, but in a much more
efficient way.

This refactoring is based on the fact that for a given neutral species n, the only variables
that are changing with each ionization collision are the cross-section of the collision σ and
the secondary_e_spectra and primary_e_spectra vectors. The rest of the calculation is the
same for all ionization collisions of the same neutral species.
Hence, we can split the ionization calculation into two parts:
- one that is the same for all ionization of the same neutral species
- one that is different for each ionization collision.
This change introduces the possibility to do some factorization and greatly reduce the
number of calculations and array memory access. Let see more in details how this works.

In this new version, we have *for each neutral species*,
1. The Ionization_fragment_1, which by analogy to the Ionization_matrix example, is equal to
    Ionization_fragment_1 = n * Ie[:, :, iE]    (no isotropic redistribution here)
2. The Ionizing_fragment_1, which by analogy to the Ionizing_matrix example, is equal to
    Ionizing_fragment_1 = n * Ie[:, :, iE]
These are calculated using the `prepare_first_ionization_fragment!()` function.

3. The Ionization_fragment_2, which is equal to
    Ionization_fragment_2 = ∑ secondary_e_spectra[i_level] .* σ[i_level]
    where the sum is done over the ionization collisions of the current neutral species
4. The Ionizing_fragment_2, which is equal to
    Ionizing_fragment_2 = ∑ primary_e_spectra[i_level] .* σ[i_level]
    where the sum is done over the ionization collisions of the current neutral species
These are calculated using the `prepare_second_ionization_fragment!()` function. This is
where the sweet factorization happens.

Then we are simply able to add the fragments to Q using the following loop (from the
function `add_ionization_fragments!()`), which needs to be run only once for the current energy.
```julia
for iI in 1:(iE - 1)
    for j in axes(Q, 2)
        for k in axes(Q, 1)
            Q[k, j, iI] += Ionization_fragment_1[1][k, j] * Ionization_fragment_2[1][iI] +  # secondaries of species 1
                            Ionization_fragment_1[2][k, j] * Ionization_fragment_2[2][iI] +  # secondaries of species 2
                            Ionization_fragment_1[3][k, j] * Ionization_fragment_2[3][iI] +  # secondaries of species 3
                            Ionizing_fragment_1[1][k, j] * Ionizing_fragment_2[1][iI] +  # primaries of species 1
                            Ionizing_fragment_1[2][k, j] * Ionizing_fragment_2[2][iI] +  # primaries of species 2
                            Ionizing_fragment_1[3][k, j] * Ionizing_fragment_2[3][iI]    # primaries of species 3
        end
    end
end
```
=#
function prepare_first_ionization_fragment!(Ionization_fragment_1, Ionizing_fragment_1,
                                            n, Ie, t, h_atm, μ_center, BeamWeight, iE)
    n_repeated_over_μt = repeat(n, length(μ_center), length(t))
    n_repeated_over_t = repeat(n, 1, length(t))

    # PRIMARY ELECTRONS
    Ionizing_fragment_1 .= n_repeated_over_μt .* @view(Ie[:, :, iE]);

    # SECONDARY ELECTRONS (ISOTROPIC)
    @views for i_μ1 in eachindex(μ_center)
        for i_μ2 in eachindex(μ_center)
            Ionization_fragment_1[(i_μ1 - 1) * length(h_atm) .+ (1:length(h_atm)), :] .+=
                n_repeated_over_t .*
                    (Ie[(i_μ2 - 1) * length(h_atm) .+ (1:length(h_atm)), :, iE]) .*
                    BeamWeight[i_μ1] ./ sum(BeamWeight)
        end
    end
end

function prepare_second_ionization_fragment!(Ionization_fragment_2, Ionizing_fragment_2,
                                             σ, E_levels, cascading, E, dE, iE)
    # Loop through the different collisions for the current neutral species
    for i_level in axes(E_levels, 1)[2:end]
        # Continue with the ionizing collisions (will produce secondary e-)
        if E_levels[i_level, 2] > 0
            # Find the energy bins where the e- in the current energy bin will degrade when
            # losing E_levels[i_level, 1] eV
            i_degrade = intersect(findall(x -> x > E[iE] - E_levels[i_level, 1], E + dE),     # find lowest bin
                                  findall(x -> x < E[iE] + dE[iE] - E_levels[i_level,1], E))  # find highest bin

            if !isempty(i_degrade) && i_degrade[1] < iE
                # Calculate the spectra of the secondary e-
                secondary_e_spectra = cascading(E, E[iE], E_levels[i_level, 1], "s");
                # We use the average energy of the e- in the current energy bin
                secondary_e_spectra = (secondary_e_spectra .+ secondary_e_spectra[[2:end; end]]) .* dE / 2

                # Calculate the distribution of the ionizing (= primary) e-, that have lost the
                # corresponding amount of energy
                primary_e_spectra = cascading(E, E[iE], E_levels[i_level, 1], "c");

                if sum(secondary_e_spectra) > 0
                    # normalise
                    secondary_e_spectra = secondary_e_spectra ./ sum(secondary_e_spectra)
                    # and multiply with the number of e- created
                    secondary_e_spectra = E_levels[i_level, 2] .* secondary_e_spectra
                end
                if sum(primary_e_spectra) > 0
                    # normalise
                    primary_e_spectra = primary_e_spectra ./ sum(primary_e_spectra)
                end

                Ionization_fragment_2 .+= secondary_e_spectra .* σ[i_level, iE]
                Ionizing_fragment_2 .+= primary_e_spectra .* σ[i_level, iE]
            end
        end
    end
end

function add_ionization_fragments!(Q, iE,
                                   Ionization_fragment_1, Ionizing_fragment_1,
                                   Ionization_fragment_2, Ionizing_fragment_2)
    @tturbo for iI in 1:(iE - 1)
        for j in axes(Q, 2)
            for k in axes(Q, 1)
                Q[k, j, iI] += Ionization_fragment_1[1][k, j] * Ionization_fragment_2[1][iI] +  # secondaries
                               Ionization_fragment_1[2][k, j] * Ionization_fragment_2[2][iI] +  # secondaries
                               Ionization_fragment_1[3][k, j] * Ionization_fragment_2[3][iI] +  # secondaries
                               Ionizing_fragment_1[1][k, j] * Ionizing_fragment_2[1][iI] +  # primaries
                               Ionizing_fragment_1[2][k, j] * Ionizing_fragment_2[2][iI] +  # primaries
                               Ionizing_fragment_1[3][k, j] * Ionizing_fragment_2[3][iI]    # primaries
            end
        end
    end
end
