using MAT
using ProgressMeter
using Dates


function calculate_e_transport(altitude_max, θ_lims, E_max, B_angle_to_zenith, t, n_loop, input_file)

    # Get atmosphere
    println("Calling Matlab for the setup...")
    h_atm, n_neutrals, ne, Te, E, dE, 
        E_levels_neutrals, σ_neutrals, secondary_e,
        θ_lims, μ_lims, μ_center, μ_scatterings = setup(altitude_max, θ_lims, E_max);

    # Initialise
    Q  = zeros(length(h_atm) * length(μ_center), length(t), length(E));
    Ie = zeros(length(h_atm) * length(μ_center), length(t), length(E));
    I0 = zeros(length(h_atm) * length(μ_center), length(E));    # starting e- flux profile

    # Load incoming flux
    Ie_top = Ie_top_from_file(input_file, μ_center, t, E, n_loop);

    # Make a finer θ for the scattering calculations
    finer_θ = range(0, π, length=721);
    # Calculate the phase functions and put them in a Tuple
    phaseN2e, phaseN2i = phase_fcn_N2(finer_θ, E);
    phaseO2e, phaseO2i = phase_fcn_O2(finer_θ, E);
    phaseOe, phaseOi = phase_fcn_O(finer_θ, E);
    phase_fcn_neutrals = ((phaseN2e, phaseN2i), (phaseO2e, phaseO2i), (phaseOe, phaseOi));
    cascading_neutrals = (cascading_N2, cascading_O2, cascading_O)

    # Create the folder to save the data to
    savedir = string(pkgdir(Aurora, "data"), "/",
                            Dates.format(now(), "yyyymmdd-HHMM"))
    println("")
    println("Results will be save at ", savedir)
    mkdir(savedir)
    # And save the simulation parameters in it
    save_parameters(altitude_max, θ_lims, E_max, B_angle_to_zenith, t, n_loop, input_file, savedir)




    ## Looping over n_loop
    for i in 1:n_loop

        D = make_D(E, dE, θ_lims);
        # Extract the top flux for the current loop
        Ie_top_local = Ie_top[:, (1 + (i - 1) * (length(t) - 1)) : (length(t) + (i - 1) * (length(t) - 1)), :]

        p = Progress(length(E), desc=string("Calculating flux for loop ", i, "/", n_loop))
        # Looping over energy
        for iE in length(E):-1:1
            A = make_A(n_neutrals, σ_neutrals, ne, Te, E, dE, iE);
            
            B, B2B_inelastic_neutrals = make_B(n_neutrals, σ_neutrals, E_levels_neutrals, 
                                                phase_fcn_neutrals, dE, iE, μ_scatterings.Pmu2mup,
                                                μ_scatterings.BeamWeight_relative, finer_θ);

            # Compute the flux of e-
            Ie[:, :, iE] = Crank_Nicolson_Optimized(t, h_atm ./ cosd(B_angle_to_zenith), μ_center, v_of_E(E[iE]), 
                                            A, B, D[iE, :], Q[:, :, iE], Ie_top_local[:, :, iE], I0[:, iE]);

            # Update the cascading of e-
            update_Q!(Q, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals, B2B_inelastic_neutrals, 
                        cascading_neutrals, E, dE, iE, μ_scatterings.BeamWeight_discrete, μ_center)

            next!(p)
        end

        # Update the starting e- flux profile
        I0 = Ie[:, end, :]

        # Save results for the n_loop
        save_results(Ie, E, t, μ_lims, h_atm, I0, μ_scatterings, n_loop, savedir, i)
    end


end