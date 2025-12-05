using Dates: Dates, now
using ProgressMeter: Progress, next!
using Term: @bold

using SparseArrays: SparseArrays, spzeros, spdiagm
using KLU: KLU, klu

# Practical cache storage
mutable struct Cache
    # Fields for steady-state scheme
    KLU::KLU.KLUFactorization{Float64, Int64}
    Mlhs::SparseArrays.SparseMatrixCSC{Float64, Int64}
    mapping::Matrix{Dict{Tuple{Symbol,Int},Int}}
    Ddz_Up::SparseArrays.SparseMatrixCSC{Float64, Int64}
    Ddz_Down::SparseArrays.SparseMatrixCSC{Float64, Int64}
    Ddiffusion::SparseArrays.SparseMatrixCSC{Float64, Int64}
    # Additional fields for Crank-Nicolson scheme
    Mrhs::SparseArrays.SparseMatrixCSC{Float64, Int64}
    mapping_lhs::Matrix{Dict{Tuple{Symbol,Int},Int}}
    mapping_rhs::Matrix{Dict{Tuple{Symbol,Int},Int}}
    # Cached repeated density matrices for ionization calculations (one per species)
    n_repeated_over_μt::Vector{Matrix{Float64}}
    n_repeated_over_t::Vector{Matrix{Float64}}
    # Cached temporary matrix for scattered flux calculation (inelastic collisions)
    Ie_scatter::Matrix{Float64}
    # Cached ionization fragment arrays (one per species)
    Ionization_fragment_1::Vector{Matrix{Float64}}
    Ionizing_fragment_1::Vector{Matrix{Float64}}
    Ionization_fragment_2::Vector{Vector{Float64}}
    Ionizing_fragment_2::Vector{Vector{Float64}}
end

function Cache()
    KLU = klu(spdiagm(0 => ones(1)))
    Mlhs = spdiagm(0 => ones(1))
    Mrhs = spdiagm(0 => ones(1))
    Ddz_Up = spdiagm(0 => ones(1))
    Ddz_Down = spdiagm(0 => ones(1))
    Ddiffusion = spdiagm(0 => ones(1))
    mapping = Matrix{Dict{Tuple{Symbol,Int},Int}}(undef, 1, 1)
    mapping[1, 1] = Dict{Tuple{Symbol,Int},Int}()
    mapping_lhs = Matrix{Dict{Tuple{Symbol,Int},Int}}(undef, 1, 1)
    mapping_lhs[1, 1] = Dict{Tuple{Symbol,Int},Int}()
    mapping_rhs = Matrix{Dict{Tuple{Symbol,Int},Int}}(undef, 1, 1)
    mapping_rhs[1, 1] = Dict{Tuple{Symbol,Int},Int}()

    n_repeated_over_μt = [zeros(1, 1) for _ in 1:3]  # 3 species (N2, O2, O)
    n_repeated_over_t = [zeros(1, 1) for _ in 1:3]
    Ie_scatter = zeros(1, 1)
    Ionization_fragment_1 = [zeros(1, 1) for _ in 1:3]
    Ionizing_fragment_1 = [zeros(1, 1) for _ in 1:3]
    Ionization_fragment_2 = [zeros(1) for _ in 1:3]
    Ionizing_fragment_2 = [zeros(1) for _ in 1:3]
    return Cache(KLU, Mlhs, mapping, Ddz_Up, Ddz_Down, Ddiffusion, Mrhs, mapping_lhs, mapping_rhs,
                 n_repeated_over_μt, n_repeated_over_t, Ie_scatter,
                 Ionization_fragment_1, Ionizing_fragment_1,
                 Ionization_fragment_2, Ionizing_fragment_2)
end

# Constructor that pre-fills the density matrices
function Cache(n_neutrals, n_μ::Int, n_t::Int, n_z::Int, n_E::Int)
    KLU = klu(spdiagm(0 => ones(1)))

    number_species = length(n_neutrals)

    # Pre-allocate and fill the density matrices for all species
    n_repeated_over_μt = Vector{Matrix{Float64}}(undef, number_species)
    n_repeated_over_t = Vector{Matrix{Float64}}(undef, number_species)

    for i_species in 1:number_species
        n = n_neutrals[i_species]

        # Allocate and fill n_repeated_over_μt (size: n_z × n_μ, n_t)
        n_repeated_over_μt[i_species] = Matrix{Float64}(undef, n_z * n_μ, n_t)
        for i_t in 1:n_t
            for i_μ in 1:n_μ
                @views n_repeated_over_μt[i_species][(i_μ - 1) * n_z .+ (1:n_z), i_t] .= n
            end
        end

        # Allocate and fill n_repeated_over_t (size: n_z, n_t)
        n_repeated_over_t[i_species] = Matrix{Float64}(undef, n_z, n_t)
        for i_t in 1:n_t
            @views n_repeated_over_t[i_species][:, i_t] .= n
        end
    end

    # Pre-allocate Ie_scatter matrix with correct size (n_z × n_μ, n_t)
    Ie_scatter = Matrix{Float64}(undef, n_z * n_μ, n_t)

    # Pre-allocate ionization fragment arrays for all species
    Ionization_fragment_1 = Vector{Matrix{Float64}}(undef, number_species)
    Ionizing_fragment_1 = Vector{Matrix{Float64}}(undef, number_species)
    Ionization_fragment_2 = Vector{Vector{Float64}}(undef, number_species)
    Ionizing_fragment_2 = Vector{Vector{Float64}}(undef, number_species)

    for i_species in 1:number_species
        Ionization_fragment_1[i_species] = zeros(n_z * n_μ, n_t)
        Ionizing_fragment_1[i_species] = zeros(n_z * n_μ, n_t)
        Ionization_fragment_2[i_species] = zeros(n_E)
        Ionizing_fragment_2[i_species] = zeros(n_E)
    end

    # Create minimal placeholder matrices with concrete types
    Mlhs = spdiagm(0 => ones(1))
    Mrhs = spdiagm(0 => ones(1))
    Ddz_Up = spdiagm(0 => ones(1))
    Ddz_Down = spdiagm(0 => ones(1))
    Ddiffusion = spdiagm(0 => ones(1))
    mapping = Matrix{Dict{Tuple{Symbol,Int},Int}}(undef, 1, 1)
    mapping[1, 1] = Dict{Tuple{Symbol,Int},Int}()
    mapping_lhs = Matrix{Dict{Tuple{Symbol,Int},Int}}(undef, 1, 1)
    mapping_lhs[1, 1] = Dict{Tuple{Symbol,Int},Int}()
    mapping_rhs = Matrix{Dict{Tuple{Symbol,Int},Int}}(undef, 1, 1)
    mapping_rhs[1, 1] = Dict{Tuple{Symbol,Int},Int}()

    return Cache(KLU, Mlhs, mapping, Ddz_Up, Ddz_Down, Ddiffusion, Mrhs, mapping_lhs, mapping_rhs,
                 n_repeated_over_μt, n_repeated_over_t, Ie_scatter,
                 Ionization_fragment_1, Ionizing_fragment_1,
                 Ionization_fragment_2, Ionizing_fragment_2)
end

function calculate_e_transport(altitude_lims, θ_lims, E_max, B_angle_to_zenith, t_total, dt,
    msis_file, iri_file, savedir, INPUT_OPTIONS, CFL_number = 64;
    n_loop = nothing, max_memory_gb = 8)

    ## Get atmosphere
    h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
    μ_scatterings = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file);

    ## CFL criteria
    # The time grid over which the simulation is run needs to be fine enough to ensure that
    # the results are correct. Here we check the CFL criteria and reduce the time grid
    # accordingly
    t, CFL_factor = CFL_criteria(t_total, dt, h_atm, v_of_E(E_max), CFL_number)

    ## Calculate n_loop automatically if not provided
    if isnothing(n_loop)
        n_loop = calculate_n_loop(t, length(h_atm), length(μ_center), length(E);
                                  max_memory_gb = max_memory_gb)
    else
        # User provided n_loop, warn if it smaller than recommended
        recommended_n_loop = calculate_n_loop(t, length(h_atm), length(μ_center), length(E);
                                              max_memory_gb = max_memory_gb, verbose = true)
        if n_loop < recommended_n_loop
            @warn "Provided n_loop=$n_loop is smaller than recommended ($recommended_n_loop). " *
                  "This may cause out of memory issues and crashes."
        end
    end

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
    save_parameters(altitude_lims, θ_lims, E_max, B_angle_to_zenith, t_total, dt, t, n_loop,
        CFL_number, INPUT_OPTIONS, savedir)
    save_neutrals(h_atm, n_neutrals, ne, Te, Tn, savedir)

    ## Calculate the number of time steps per loop
    # Each loop processes (length(t) - 1) / n_loop + 1 time steps, with overlapping endpoints
    n_t_per_loop = (length(t) - 1) ÷ n_loop + 1
    t_per_loop = range(0, step = t[2] - t[1], length = n_t_per_loop)

    ## Initialize cache with pre-filled density matrices
    cache = Cache(n_neutrals, length(μ_center), n_t_per_loop, length(h_atm), length(E))

    ## Initialise
    I0 = zeros(length(h_atm) * length(μ_center), length(E));    # starting e- flux profile
    Ie = zeros(length(h_atm) * length(μ_center), n_t_per_loop, length(E))
    n_t_save = length(0:dt:t_per_loop[end])
    Ie_save = zeros(length(h_atm) * length(μ_center), n_t_save, length(E))

    ## Initialize transport matrices container
    matrices = initialize_transport_matrices(h_atm, μ_center, t_per_loop, E, dE, θ_lims)
    # Pre-compute the D matrix (energy × angle) and the diffusion operator (altitude x altitude)
    update_D!(matrices.D, E, dE, θ_lims)
    update_Ddiffusion!(matrices.Ddiffusion, h_atm ./ cosd(B_angle_to_zenith))
    matrices.Ddiffusion[1, 1] = 0

    ## Precalculate the B2B fragment
    B2B_fragment = prepare_beams2beams(μ_scatterings.BeamWeight_relative, μ_scatterings.Pmu2mup);

    ## Looping over n_loop
    for i_loop in 1:n_loop

        # Reset Ie and Q for each new loop
        fill!(Ie, 0.0)
        fill!(matrices.Q, 0.0)

        # Extract the top flux for the current loop
        Ie_top_local = Ie_top[:, (1 + (i_loop - 1) * (n_t_per_loop - 1)) : (n_t_per_loop + (i_loop - 1) * (n_t_per_loop - 1)), :]

        p = Progress(length(E); desc=string("Calculating flux for loop ", i_loop, "/", n_loop), dt=1.0)
        # Looping over energy
        for iE in length(E):-1:1
            # Update matrices A and B for current energy
            B2B_inelastic_neutrals = update_matrices!(matrices, n_neutrals, σ_neutrals, ne, Te,
                                                      E_levels_neutrals, phase_fcn_neutrals,
                                                      E, dE, iE, B2B_fragment, μ_scatterings.theta1)

            # Compute the flux of e-
            if iE == length(E)
                @views Crank_Nicolson_optimized!(Ie[:, :, iE], t_per_loop,
                                                 h_atm ./ cosd(B_angle_to_zenith), μ_center,
                                                 v_of_E(E[iE]), matrices, iE,
                                                 Ie_top_local[:, :, iE], I0[:, iE], cache,
                                                 first_iteration = true)
            else
                @views Crank_Nicolson_optimized!(Ie[:, :, iE], t_per_loop,
                                                 h_atm ./ cosd(B_angle_to_zenith), μ_center,
                                                 v_of_E(E[iE]), matrices, iE,
                                                 Ie_top_local[:, :, iE], I0[:, iE], cache)
            end

            # Update the cascading of e-
            update_Q!(matrices, Ie, h_atm, t_per_loop, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals,
                        B2B_inelastic_neutrals, cascading_neutrals, E, dE, iE, μ_scatterings.BeamWeight,
                        μ_center, cache)

            next!(p)
        end

        # Update the starting e- flux profile
        I0 .= Ie[:, end, :]

        # Subsample Ie into Ie_save buffer
        Ie_save .= @view Ie[:, 1:CFL_factor:end, :]
        # Save results for the n_loop
        save_results(Ie_save, E, t_per_loop, μ_lims, h_atm, I0, μ_scatterings, i_loop, CFL_factor, savedir)
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

    ## Initialize cache with pre-filled density matrices
    cache = Cache(n_neutrals, length(μ_center), 1, length(h_atm), length(E))

    ## Precalculate the B2B fragment
    B2B_fragment = prepare_beams2beams(μ_scatterings.BeamWeight_relative, μ_scatterings.Pmu2mup);

    ## Initialize transport matrices container
    matrices = initialize_transport_matrices(h_atm, μ_center, 1:1:1, E, dE, θ_lims)
    # Pre-compute the D matrix (energy × angle) and the diffusion operator (altitude x altitude)
    update_D!(matrices.D, E, dE, θ_lims)
    update_Ddiffusion!(matrices.Ddiffusion, h_atm ./ cosd(B_angle_to_zenith))
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
            @views steady_state_scheme_optimized!(Ie[:, 1, iE],
                                                  h_atm ./ cosd(B_angle_to_zenith),
                                                  μ_center, matrices, iE,
                                                  Ie_top_local[:, iE], cache,
                                                  first_iteration = true)
        else
            @views steady_state_scheme_optimized!(Ie[:, 1, iE],
                                                  h_atm ./ cosd(B_angle_to_zenith),
                                                  μ_center, matrices, iE,
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
