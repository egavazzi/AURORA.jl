function loss_to_thermal_electrons(E, ne, Te)
    kB = 1.380662e-23;      # Boltzmann constant [J/K]
    qₑ = 1.6021773e-19;    # elementary charge [C]
    if length(E) > 1
        sz_E = length(E);
        sz_ne = length(ne);

        ne = repeat(ne, outer=(1, sz_E));
        Te = repeat(Te, outer=(1, sz_E));
        E = repeat(E, outer=(1, sz_ne))';

        Ee = kB / qₑ * Te;

        Le = real(3.0271e-10 * ne .^ .97 .* ((E .- Ee) ./ (E .- 0.53 * Ee) .^ 2.36 ./ E .^ .44 / v_of_E(E)));

        Le[E .< Ee] .= 0;
    else
        Ee = 8.618e-5 * Te;

        Le = real(3.0271e-10 * ne .^ .97 .* ((E .- Ee) ./ (E .- 0.53 * Ee) .^ 2.36 ./ E .^ .44 / v_of_E(E)));
        Le[E .< Ee] .= 0;
    end
    return Le
end

function beams2beams(phase_fcn, Pmu2mup, BeamWeight_relative)
    B2B = zeros(size(Pmu2mup, 3),size(Pmu2mup, 3));
    for i = size(Pmu2mup, 3):-1:1
        B2B[i, :] = BeamWeight_relative * (@view(Pmu2mup[:, :, i]) * phase_fcn);
    end
    return B2B
end

## ----------------------------------------------------- ##

function make_A(n_neutrals, σ_neutrals, ne, Te, E, dE, iE)
    A = zeros(length(n_neutrals[1]));
    # Loop over the neutral species
    for i1 in 1:length(n_neutrals)
        n = n_neutrals[i1];  # Neutral density
        σ = σ_neutrals[i1];  # Array with collision cross sections

        # add elastic collisions
        A = A + n * σ[1, iE];

        # add inelastic and ionization collisions
        for i2 in 2:size(σ, 1)      # Loop over the different collisions, because
            A = A + n * σ[i2, iE];  # they have different cross sections
        end

        # add losses due to electron-electron collisions
        A = A + loss_to_thermal_electrons(E[iE], ne, Te) / dE[iE];
    end
    return A
end

function make_B(n_neutrals, σ_neutrals, E_levels_neutrals, phase_fcn_neutrals, dE, iE, Pmu2mup, BeamWeight_relative, finer_θ)
    B = zeros(length(n_neutrals[1]), size(Pmu2mup, 3), size(Pmu2mup, 3));
    # Loop over the neutral species
    for i in 1:length(n_neutrals)
        n = n_neutrals[i];                  # Neutral density
        σ = σ_neutrals[i];                  # Array with collision cross sections
        E_levels = E_levels_neutrals[i];    # Array with collision enery levels and number of secondary e-
        phase_fcn = phase_fcn_neutrals[i];  # Tuple with two phase function arrays, the first for elastic collisions 
                                            # and the second for inelastic collisions

        # Convert to 3D the scattering probabilities that are in 1D
        phase_fcn_e = convert_phase_fcn_to_3D(phase_fcn[1], finer_θ);
        phase_fcn_i = convert_phase_fcn_to_3D(phase_fcn[2], finer_θ);
        B2B_elastic = beams2beams(phase_fcn_e, Pmu2mup, BeamWeight_relative);
        B2B_inelastic = beams2beams(phase_fcn_i, Pmu2mup, BeamWeight_relative);

        # add scattering from elastic collisions
        for i1 in size(B2B_elastic, 1):-1:1
            for i2 in size(B2B_elastic, 2):-1:1
                B[:, i1, i2] .= @view(B[:, i1, i2]) .+ n .* σ[1, iE] .* B2B_elastic[i1, i2];
            end
        end

        # add scattering from inelastic and ionization collisions
        for i1 in 2:size(σ, 1)
            for i2 in size(B2B_inelastic, 1):-1:1
                for i3 in size(B2B_inelastic, 2):-1:1
                    # The second factor corrects for the case where the energy loss
                    # E_levels[i1,1] is smaller than the width in energy in the energy bin.
                    # That is, when dE[iE] > E_levels[i1,1], only the fraction
                    # E_levels[i1,1]/dE is lost from the energy bin [E[iE], E[iE] + dE[iE]].
                    B[:, i2, i3] .= @view(B[:, i2, i3]) .+ n .* σ[i1, iE] .* B2B_inelastic[i2, i3] .* 
                                                    max(0, 1 - E_levels[i1, 1] ./ dE[iE]);
                end
            end
        end
    end
    return B
end

function make_D(E, dE, θ_lims)
    θ_lims = deg2rad.(θ_lims)
    nE = 3
    nθ = 3
    # n_ti = 701
    # n_thi = 401
    D_e = zeros(length(E), length(θ_lims) -1)
    for iE in length(E):-1:1
        v = range(v_of_E(E[iE]), v_of_E(E[iE] + dE[iE]), length=nE)
        for iθ in 1:(length(θ_lims) - 1)
            θa = θ_lims[iθ]
            θb = θ_lims[iθ + 1]
            if θ_lims[iθ] == π/2
                θa = θ_lims[iθ] * 0.8 + 0.2 * θ_lims[iθ + 1]
            end
            if θ_lims[iθ + 1] == π/2
                θb = θ_lims[iθ] * 0.2 + 0.8 * θ_lims[iθ + 1]
            end
            θ = range(θa, θb, length=nθ)
            # θ4i = range(minimum(θ), maximum(θ), n_thi)
            v_par = [A * cos(B) for A in v, B in θ]
            t_arrival = 500e3 ./ v_par
            at_a = (maximum(t_arrival) + minimum(t_arrival)) / 2
            dt_a = (maximum(t_arrival) - minimum(t_arrival))
            D = (dt_a / 4)^2 / at_a
            
            D_e[iE, iθ] = abs(D)
        end
    end
    return D_e
end