using SparseArrays

function make_big_B2B_matrix(B2B_inelastic, n, h_atm)
    idx1 = Vector{Float64}(undef, 0)
    idx2 = Vector{Float64}(undef, 0)
    aB2B = Vector{Float64}(undef, 0)
    for i1 in axes(B2B_inelastic, 1)
        for i2 in axes(B2B_inelastic, 2)
            append!(idx1, (i1 - 1) * length(h_atm) .+ (1:length(h_atm)))
            append!(idx2, (i2 - 1) * length(h_atm) .+ (1:length(h_atm)))
            append!(aB2B, n * B2B_inelastic[i1, i2])
        end
    end
    AB2B = sparse(idx1, idx2, aB2B)
    return AB2B
end

function add_inelastic_collisions!(Q, Ie, h_atm, n, σ, E_levels, B2B_inelastic, E, dE, iE)
    Ie_degraded = Matrix{Float64}(undef, size(Ie, 1), size(Ie, 2))

    # Multiply each element of B2B with n (density vector) and resize to get a matrix that can
    # be multiplied with Ie
    AB2B =  make_big_B2B_matrix(B2B_inelastic, n, h_atm)
    Ie_scatter = AB2B * @view(Ie[:, :, iE])

    # Loop over the energy levels of the collisions with the i-th neutral specie
    for i_level in 2:size(E_levels, 1)
        if E_levels[i_level, 2] <= 0  # these collisions should not produce secondary e-
            # The flux of e- degraded from energy bin [E[iE], E[iE] + dE[iE]] to any lower energy
            # bin by excitation of the E_levels[i_level] state of the current specie.
            # The second factor corrects for the case where the energy loss is maller than the width
            # in energy of the energy bin. That is, when dE[iE] > E_levels[i_level,1], only the
            # fraction E_levels[i_level,1] / dE[iE] is lost from the energy bin [E[iE], E[iE] + dE[iE]].

            Ie_degraded .= (σ[i_level, iE] * min(1, E_levels[i_level, 1] / dE[iE])) .* Ie_scatter

            # Find the energy bins where the e- in the current energy bin will degrade when losing
            # E_levels[i_level, 1] eV
            i_degrade = intersect(  findall(x -> x > E[iE] - E_levels[i_level, 1], E + dE),     # find lowest bin
                                    findall(x -> x < E[iE] + dE[iE] - E_levels[i_level,1], E))  # find highest bin

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
                Threads.@threads for i_u in findall(x -> x != 0, partition_fraction)
                    @view(Q[:, :, i_degrade[i_u]]) .+=  max.(0, Ie_degraded) .* partition_fraction[i_u]
                end
            end
        end
    end
end

using Polyester
function add_ionization_collisions!(Q, Ie, h_atm, t, n, σ, E_levels, cascading, E, dE, iE, BeamWeight, μ_center)
    Ionization = Matrix{Float64}(undef, size(Ie, 1), size(Ie, 2))
    Ionizing = Matrix{Float64}(undef, size(Ie, 1), size(Ie, 2))
    n_repeated_over_μt = repeat(n, length(μ_center), length(t))
    n_repeated_over_t = repeat(n, 1, length(t))

    for i_level = 2:size(E_levels, 1)
        if E_levels[i_level, 2] > 0    # these collisions should produce secondary e-
            # Find the energy bins where the e- in the current energy bin will degrade when losing
            # E_levels[i_level, 1] eV
            i_degrade = intersect(findall(x -> x > E[iE] - E_levels[i_level, 1], E + dE),     # find lowest bin
                                  findall(x -> x < E[iE] + dE[iE] - E_levels[i_level,1], E))  # find highest bin

            if !isempty(i_degrade) && i_degrade[1] < iE
                # ISOTROPIC SECONDARY ELECTRONS

                Ionizing .= n_repeated_over_μt .* (σ[i_level, iE] .* @view(Ie[:, :, iE]));
                @views for i_μ in eachindex(μ_center)
                    Ionization[(i_μ - 1) * length(h_atm) .+ (1:length(h_atm)), :] .=
                        max.(0, n_repeated_over_t .*
                             (σ[i_level, iE] .*
                             Ie[(i_μ - 1) * length(h_atm) .+ (1:length(h_atm)), :, iE]) .*
                             BeamWeight[i_μ] ./ sum(BeamWeight))
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

                # This is some kind of check that it is not empty?
                e_ionized_distribution = max.(secondary_e_spectra, primary_e_spectra)
                e_ionized_distribution[isnan.(e_ionized_distribution)] .= 0


                # and finally add this to the flux of degrading e-
                Nthreads = 6 # will run on 6 threads, which seems to be optimal on my machine
                nbatch = Int(floor(iE / Nthreads)) - 1
                @batch minbatch=nbatch for iI in 1:(iE - 1)
                    @view(Q[:, :, iI]) .+= Ionization .* secondary_e_spectra[iI] .+
                                           Ionizing .* primary_e_spectra[iI]
                end
            end
        end
    end
end

function update_Q!(Q, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals, B2B_inelastic_neutrals,
                    cascading_neutrals, E, dE, iE, BeamWeight, μ_center)

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
        cascading = cascading_neutrals[i];          # Cascading function for the current i-th specie

        add_inelastic_collisions!(Q, Ie, h_atm, n, σ, E_levels, B2B_inelastic, E, dE, iE)
        add_ionization_collisions!(Q, Ie, h_atm, t, n, σ, E_levels, cascading, E, dE, iE, BeamWeight, μ_center)
    end
end
