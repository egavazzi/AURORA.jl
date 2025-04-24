using Dates: Dates, now
using ProgressMeter: Progress, next!
using Term: @bold

function calculate_e_transport(altitude_lims, θ_lims, E_max, B_angle_to_zenith, t_sampling,
    n_loop, msis_file, iri_file, savedir, INPUT_OPTIONS, CFL_number = 64)

    ## Get atmosphere
    h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
    μ_scatterings = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file);

    ## Initialise
    I0 = zeros(length(h_atm) * length(μ_center), length(E));    # starting e- flux profile

    ## CFL criteria
    # The time grid over which the simulation is run needs to be fine enough to ensure that
    # the results are correct. Here we check the CFL criteria and reduce the time grid
    # accordingly
    t, CFL_factor = CFL_criteria(t_sampling, h_atm, v_of_E(E_max), CFL_number)
    # TODO: fix properly the time sampling of incoming data from file

    ## Load incoming flux
    if INPUT_OPTIONS.input_type == "from_old_matlab_file"
        Ie_top = Ie_top_from_old_matlab_file(t, E, n_loop, μ_center, INPUT_OPTIONS.input_file);
    elseif INPUT_OPTIONS.input_type == "from_file"
        Ie_top = Ie_top_from_file(t, E, μ_center, n_loop, INPUT_OPTIONS.input_file)
    elseif INPUT_OPTIONS.input_type == "flickering"
        Ie_top = Ie_top_flickering(t, E, dE, n_loop, μ_center, h_atm,
                                    μ_scatterings.BeamWeight, INPUT_OPTIONS.IeE_tot,
                                    INPUT_OPTIONS.z₀, INPUT_OPTIONS.E_min, INPUT_OPTIONS.f,
                                    INPUT_OPTIONS.Beams, INPUT_OPTIONS.modulation)
    elseif INPUT_OPTIONS.input_type == "constant_onset"
        Ie_top = Ie_top_constant(t, E, dE, n_loop, μ_center, h_atm,
                                μ_scatterings.BeamWeight, INPUT_OPTIONS.IeE_tot,
                                INPUT_OPTIONS.z₀, INPUT_OPTIONS.E_min, INPUT_OPTIONS.Beams,
                                INPUT_OPTIONS.t0, INPUT_OPTIONS.t1)
    elseif INPUT_OPTIONS.input_type == "from_ketchup_file"
        Ie_top = Ie_top_from_ketchup(t, E, n_loop, μ_center, INPUT_OPTIONS.input_file);
    end

    ## Calculate the phase functions and put them in a Tuple
    phaseN2e, phaseN2i = phase_fcn_N2(μ_scatterings.theta1, E);
    phaseO2e, phaseO2i = phase_fcn_O2(μ_scatterings.theta1, E);
    phaseOe, phaseOi = phase_fcn_O(μ_scatterings.theta1, E);
    phase_fcn_neutrals = ((phaseN2e, phaseN2i), (phaseO2e, phaseO2i), (phaseOe, phaseOi));
    cascading_neutrals = (cascading_N2, cascading_O2, cascading_O) # tuple of functions

    ## And save the simulation parameters in it
    save_parameters(altitude_lims, θ_lims, E_max, B_angle_to_zenith, t_sampling, t, n_loop,
        CFL_number, INPUT_OPTIONS, savedir)
    save_neutrals(h_atm, n_neutrals, ne, Te, Tn, savedir)

    # Initialize arrays for the ionization collisions part of the energy degradation
    Ionization_matrix = [zeros(length(h_atm) * length(μ_center), length(t)) for _ in 1:15]
    Ionizing_matrix = [zeros(length(h_atm) * length(μ_center), length(t)) for _ in 1:15]
    secondary_vector = [zeros(length(E)) for _ in 1:15]
    primary_vector = [zeros(length(E)) for _ in 1:15]

    ## Looping over n_loop
    for i in 1:n_loop
        Q  = zeros(length(h_atm) * length(μ_center), length(t), length(E));
        Ie = zeros(length(h_atm) * length(μ_center), length(t), length(E));

        D = make_D(E, dE, θ_lims);
        # Extract the top flux for the current loop
        Ie_top_local = Ie_top[:, (1 + (i - 1) * (length(t) - 1)) : (length(t) + (i - 1) * (length(t) - 1)), :]

        p = Progress(length(E); desc=string("Calculating flux for loop ", i, "/", n_loop), dt=1.0)
        # Looping over energy
        for iE in length(E):-1:1
            A = make_A(n_neutrals, σ_neutrals, ne, Te, E, dE, iE);

            B, B2B_inelastic_neutrals = make_B(n_neutrals, σ_neutrals, E_levels_neutrals,
                                                phase_fcn_neutrals, dE, iE, μ_scatterings.Pmu2mup,
                                                μ_scatterings.BeamWeight_relative, μ_scatterings.theta1);

            # Compute the flux of e-
            Ie[:, :, iE] = Crank_Nicolson(t, h_atm ./ cosd(B_angle_to_zenith), μ_center, v_of_E(E[iE]),
                                            A, B, D[iE, :], Q[:, :, iE], Ie_top_local[:, :, iE], I0[:, iE]);

            # Update the cascading of e-
            update_Q_turbo!(Q, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals,
                        B2B_inelastic_neutrals, cascading_neutrals, E, dE, iE,
                        μ_scatterings.BeamWeight, μ_center,
                        Ionization_matrix, Ionizing_matrix, secondary_vector, primary_vector)

            next!(p)
        end

        # Update the starting e- flux profile
        I0 = Ie[:, end, :]

        # Save results for the n_loop
        save_results(Ie, E, t, μ_lims, h_atm, I0, μ_scatterings, i, CFL_factor, savedir)
    end

end





function calculate_e_transport_steady_state(altitude_lims, θ_lims, E_max, B_angle_to_zenith,
    msis_file, iri_file, savedir, INPUT_OPTIONS)
    ## Get atmosphere
    h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
    μ_scatterings = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file);

    ## Initialise
    I0 = zeros(length(h_atm) * length(μ_center), length(E)); # starting e- flux profile

    ## Load incoming flux
    if INPUT_OPTIONS.input_type == "from_old_matlab_file"
        Ie_top = Ie_top_from_old_matlab_file(1:1:1, E, 1, μ_center, INPUT_OPTIONS.input_file);
    elseif INPUT_OPTIONS.input_type == "from_file"
        Ie_top = Ie_top_from_file(1:1:1, E, μ_center, 1, INPUT_OPTIONS.input_file)
    elseif INPUT_OPTIONS.input_type == "flickering"
        Ie_top = Ie_top_flickering(1:1:1, E, dE, 1, μ_center, h_atm,
                                    μ_scatterings.BeamWeight, INPUT_OPTIONS.IeE_tot,
                                    INPUT_OPTIONS.z₀, INPUT_OPTIONS.E_min, INPUT_OPTIONS.f,
                                    INPUT_OPTIONS.Beams, INPUT_OPTIONS.modulation)
    elseif INPUT_OPTIONS.input_type == "constant_onset"
        Ie_top = Ie_top_constant(1:1:1, E, dE, 1, μ_center, h_atm,
                                μ_scatterings.BeamWeight, INPUT_OPTIONS.IeE_tot,
                                INPUT_OPTIONS.z₀, INPUT_OPTIONS.E_min, INPUT_OPTIONS.Beams,
                                INPUT_OPTIONS.t0, INPUT_OPTIONS.t1)
    elseif INPUT_OPTIONS.input_type == "LET"
        Ie_top = Ie_with_LET(INPUT_OPTIONS.E0, INPUT_OPTIONS.Q, E, μ_center,
                             INPUT_OPTIONS.Beams, INPUT_OPTIONS.low_energy_tail)
    elseif INPUT_OPTIONS.input_type == "from_ketchup_file"
        Ie_top = Ie_top_from_ketchup(1:1:1, E, 1, μ_center, INPUT_OPTIONS.input_file);
    elseif INPUT_OPTIONS.input_type == "bgu_custom"
        Ie_top = INPUT_OPTIONS.Ie_top
    end

    ## Calculate the phase functions and put them in a Tuple
    phaseN2e, phaseN2i = phase_fcn_N2(μ_scatterings.theta1, E);
    phaseO2e, phaseO2i = phase_fcn_O2(μ_scatterings.theta1, E);
    phaseOe, phaseOi = phase_fcn_O(μ_scatterings.theta1, E);
    phase_fcn_neutrals = ((phaseN2e, phaseN2i), (phaseO2e, phaseO2i), (phaseOe, phaseOi));
    cascading_neutrals = (cascading_N2, cascading_O2, cascading_O) # tuple of functions

    ## And save the simulation parameters in it
    save_parameters(altitude_lims, θ_lims, E_max, B_angle_to_zenith, 1:1:1, 1:1:1, 1,
        0, INPUT_OPTIONS, savedir)
    save_neutrals(h_atm, n_neutrals, ne, Te, Tn, savedir)

    # Initialize arrays for the ionization collisions part of the energy degradation
    Ionization_matrix = [zeros(length(h_atm) * length(μ_center), 1) for _ in 1:15]
    Ionizing_matrix = [zeros(length(h_atm) * length(μ_center), 1) for _ in 1:15]
    secondary_vector = [zeros(length(E)) for _ in 1:15]
    primary_vector = [zeros(length(E)) for _ in 1:15]


    ## No n_loop here
    Q  = zeros(length(h_atm) * length(μ_center), 1, length(E));
    Ie = zeros(length(h_atm) * length(μ_center), 1, length(E));

    D = make_D(E, dE, θ_lims);
    # Extract the top flux for the current loop
    Ie_top_local = Ie_top[:, 1, :];

    p = Progress(length(E), desc=string("Calculating flux"))
    # Looping over energy
    for iE in length(E):-1:1
        A = make_A(n_neutrals, σ_neutrals, ne, Te, E, dE, iE);

        B, B2B_inelastic_neutrals = make_B(n_neutrals, σ_neutrals, E_levels_neutrals,
                                phase_fcn_neutrals, dE, iE, μ_scatterings.Pmu2mup,
                                μ_scatterings.BeamWeight_relative, μ_scatterings.theta1);

        # Compute the flux of e-
        Ie[:, 1, iE] = steady_state_scheme(h_atm ./ cosd(B_angle_to_zenith), μ_center, A, B,
                                    D[iE, :],  Q[:, 1, iE], Ie_top_local[:, iE])

        # Update the cascading of e-
        update_Q_turbo!(Q, Ie, h_atm, 1:1:1, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals, B2B_inelastic_neutrals,
                    cascading_neutrals, E, dE, iE, μ_scatterings.BeamWeight, μ_center,
                    Ionization_matrix, Ionizing_matrix, secondary_vector, primary_vector)

        next!(p)
    end

    save_results(Ie, E, 1:1:1, μ_lims, h_atm, I0, μ_scatterings, 1, 1, savedir)

end
