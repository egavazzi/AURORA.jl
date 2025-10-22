using Dates: Dates, now
using ProgressMeter: Progress, next!
using Term: @bold

using SparseArrays: SparseArrays, spzeros, spdiagm
using KLU: KLU, klu

# Practical cache storage
mutable struct Cache
    AB2B::SparseArrays.SparseMatrixCSC{Float64, Int64}
    KLU::KLU.KLUFactorization{Float64, Int64}
end

function Cache()
    AB2B = spzeros(1, 1)
    KLU = klu(spdiagm(0 => ones(1)))
    return Cache(AB2B, KLU)
end

function calculate_e_transport(altitude_lims, θ_lims, E_max, B_angle_to_zenith, t_sampling,
    n_loop, msis_file, iri_file, savedir, INPUT_OPTIONS, CFL_number = 64)

    ## Get atmosphere
    h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
    μ_scatterings = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file);

    ## CFL criteria
    # The time grid over which the simulation is run needs to be fine enough to ensure that
    # the results are correct. Here we check the CFL criteria and reduce the time grid
    # accordingly
    t, CFL_factor = CFL_criteria(t_sampling, h_atm, v_of_E(E_max), CFL_number)
    # TODO: fix properly the time sampling of incoming data from file   # What did I mean here?? /EG20250504

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
    print("Calculating the phase functions...")
    phaseN2e, phaseN2i = phase_fcn_N2(μ_scatterings.theta1, E);
    phaseO2e, phaseO2i = phase_fcn_O2(μ_scatterings.theta1, E);
    phaseOe, phaseOi = phase_fcn_O(μ_scatterings.theta1, E);
    phase_fcn_neutrals = ((phaseN2e, phaseN2i), (phaseO2e, phaseO2i), (phaseOe, phaseOi));
    cascading_neutrals = (cascading_N2, cascading_O2, cascading_O) # tuple of functions
    println(" done ✅")

    ## And save the simulation parameters in it
    save_parameters(altitude_lims, θ_lims, E_max, B_angle_to_zenith, t_sampling, t, n_loop,
        CFL_number, INPUT_OPTIONS, savedir)
    save_neutrals(h_atm, n_neutrals, ne, Te, Tn, savedir)

    ## Initialize cache
    cache = Cache()

    ## Initialise
    I0 = zeros(length(h_atm) * length(μ_center), length(E));    # starting e- flux profile

    ## Initialize transport matrices container
    matrices = initialize_transport_matrices(h_atm, μ_center, t, E, dE, θ_lims)
    # Pre-compute the D matrix (energy × angle)
    matrices.D .= make_D(E, dE, θ_lims)
    # Pre-compute the diffusion operator (altitude x altitude)
    matrices.Ddiffusion = d2M(h_atm ./ cosd(B_angle_to_zenith))
    matrices.Ddiffusion[1, 1] = 0

    ## Precalculate the B2B fragment
    B2B_fragment = prepare_beams2beams(μ_scatterings.BeamWeight_relative, μ_scatterings.Pmu2mup);

    ## Looping over n_loop
    for i in 1:n_loop
        # Initialize the ionospheric flux for the current loop
        Ie = zeros(length(h_atm) * length(μ_center), length(t), length(E));

        # Reset Q for each new loop
        fill!(matrices.Q, 0.0)

        # Extract the top flux for the current loop
        Ie_top_local = Ie_top[:, (1 + (i - 1) * (length(t) - 1)) : (length(t) + (i - 1) * (length(t) - 1)), :]

        p = Progress(length(E); desc=string("Calculating flux for loop ", i, "/", n_loop), dt=1.0)
        # Looping over energy
        for iE in length(E):-1:1
            # Update matrices A and B for current energy
            B2B_inelastic_neutrals = update_matrices!(matrices, n_neutrals, σ_neutrals, ne, Te,
                                                      E_levels_neutrals, phase_fcn_neutrals,
                                                      E, dE, iE, B2B_fragment, μ_scatterings.theta1)

            # Compute the flux of e-
            if iE == length(E)
                Ie[:, :, iE] = Crank_Nicolson(t, h_atm ./ cosd(B_angle_to_zenith), μ_center, v_of_E(E[iE]),
                                                         matrices, iE,
                                                         Ie_top_local[:, :, iE], I0[:, iE],
                                                         cache, first_iteration = true)
            else
                Ie[:, :, iE] = Crank_Nicolson(t, h_atm ./ cosd(B_angle_to_zenith), μ_center, v_of_E(E[iE]),
                                                         matrices, iE,
                                                         Ie_top_local[:, :, iE], I0[:, iE],
                                                         cache)
            end

            # Update the cascading of e-
            update_Q!(matrices, Ie, h_atm, t, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals,
                        B2B_inelastic_neutrals, cascading_neutrals, E, dE, iE, μ_scatterings.BeamWeight,
                        μ_center, cache)

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
        Ie_top = Ie_with_LET(INPUT_OPTIONS.E0, INPUT_OPTIONS.Q, E, dE, μ_center,
                             μ_scatterings.BeamWeight, INPUT_OPTIONS.Beams;
                             low_energy_tail = INPUT_OPTIONS.low_energy_tail)
    elseif INPUT_OPTIONS.input_type == "from_ketchup_file"
        Ie_top = Ie_top_from_ketchup(1:1:1, E, 1, μ_center, INPUT_OPTIONS.input_file);
    elseif INPUT_OPTIONS.input_type == "bgu_custom"
        Ie_top = INPUT_OPTIONS.Ie_top
    end

    ## Calculate the phase functions and put them in a Tuple
    print("Calculating the phase functions...")
    phaseN2e, phaseN2i = phase_fcn_N2(μ_scatterings.theta1, E);
    phaseO2e, phaseO2i = phase_fcn_O2(μ_scatterings.theta1, E);
    phaseOe, phaseOi = phase_fcn_O(μ_scatterings.theta1, E);
    phase_fcn_neutrals = ((phaseN2e, phaseN2i), (phaseO2e, phaseO2i), (phaseOe, phaseOi));
    cascading_neutrals = (cascading_N2, cascading_O2, cascading_O) # tuple of functions
    println(" done ✅")

    ## And save the simulation parameters in it
    save_parameters(altitude_lims, θ_lims, E_max, B_angle_to_zenith, 1:1:1, 1:1:1, 1,
        0, INPUT_OPTIONS, savedir)
    save_neutrals(h_atm, n_neutrals, ne, Te, Tn, savedir)

    ## Initialize cache
    cache = Cache()

    ## Precalculate the B2B fragment
    B2B_fragment = prepare_beams2beams(μ_scatterings.BeamWeight_relative, μ_scatterings.Pmu2mup);

    ## Initialize transport matrices container for steady state (1 time step)
    matrices = initialize_transport_matrices(h_atm, μ_center, 1:1:1, E, dE, θ_lims)

    matrices.D .= make_D(E, dE, θ_lims)
    # Pre-compute the diffusion operator (altitude x altitude)
    matrices.Ddiffusion = d2M(h_atm ./ cosd(B_angle_to_zenith))
    matrices.Ddiffusion[1, 1] = 0

    ## Initialize the ionospheric flux
    Ie = zeros(length(h_atm) * length(μ_center), 1, length(E));

    # Extract the top flux for the current loop
    Ie_top_local = Ie_top[:, 1, :];

    p = Progress(length(E), desc=string("Calculating flux"))
    # Looping over energy
    for iE in length(E):-1:1
        # Update matrices A and B for current energy
        B2B_inelastic_neutrals = update_matrices!(matrices, n_neutrals, σ_neutrals, ne, Te,
                                                  E_levels_neutrals, phase_fcn_neutrals,
                                                  E, dE, iE, B2B_fragment, μ_scatterings.theta1)

        # Compute the flux of e-
        if iE == length(E)
            Ie[:, 1, iE] = steady_state_scheme(h_atm ./ cosd(B_angle_to_zenith),
                                                          μ_center, matrices, iE,
                                                          Ie_top_local[:, iE],
                                                          cache, first_iteration = true)
        else
            Ie[:, 1, iE] = steady_state_scheme(h_atm ./ cosd(B_angle_to_zenith), μ_center,
                                               matrices, iE,
                                               Ie_top_local[:, iE], cache)
        end

        # Update the cascading of e-
        update_Q!(matrices, Ie, h_atm, 1:1:1, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals,
                  B2B_inelastic_neutrals, cascading_neutrals, E, dE, iE, μ_scatterings.BeamWeight,
                  μ_center, cache)

        next!(p)
    end

    save_results(Ie, E, 1:1:1, μ_lims, h_atm, I0, μ_scatterings, 1, 1, savedir)

end
