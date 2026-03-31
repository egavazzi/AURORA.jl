using ProgressMeter: Progress, next!

function require_cache(sim::AuroraSimulation)
    cache = sim.cache
    cache === nothing && error("Simulation not initialized. Call initialize!(sim) or run!(sim).")
    return cache
end


# Compute phase functions for electron and ion scattering off neutral molecules (N2, O2, O).
# These phase functions describe the angular distribution of particles after collisions.
function compute_phase_functions(model::AuroraModel)
    E_centers = model.energy_grid.E_centers
    θ_scatter = model.scattering.θ_scatter

    print("Calculating the phase functions...")
    phaseN2e, phaseN2i = phase_fcn_N2(θ_scatter, E_centers)
    phaseO2e, phaseO2i = phase_fcn_O2(θ_scatter, E_centers)
    phaseOe, phaseOi = phase_fcn_O(θ_scatter, E_centers)
    println(" done ✅")

    phase_fcn_neutrals = ((phaseN2e, phaseN2i), (phaseO2e, phaseO2i), (phaseOe, phaseOi))
    return phase_fcn_neutrals
end

function build_transport_cache(sim::AuroraSimulation)
    # Extract model geometry and grids
    model = sim.model
    z = model.altitude_grid.h
    μ_center = model.pitch_angle_grid.μ_center
    neutral_densities = n_neutrals(model.ionosphere)
    n_E = model.energy_grid.n

    # Set up time grid (single step for steady-state, multiple steps for time-dependent)
    n_t = isnothing(sim.time) ? 1 : sim.time.n_t_per_loop
    t_loop = isnothing(sim.time) ? (1:1:1) : range(0.0, step=sim.time.dt_resolved, length=sim.time.n_t_per_loop)

    # Initialize solver and physical process caches
    solver = SolverCache()
    degradation = DegradationCache(neutral_densities, length(μ_center), n_t, length(z), n_E)
    matrices = initialize_transport_matrices(model, t_loop)
    update_D!(matrices.D, model)
    update_Ddiffusion!(matrices.Ddiffusion, model)

    # Pre-compute scattering
    phase_fcn_neutrals = compute_phase_functions(model)
    B2B_fragment = prepare_beams2beams(model.scattering.Ω_subbeam_relative, model.scattering.P_scatter)

    # Pre-compute input flux data
    Ie_top = isnothing(sim.time) ? compute_flux(sim.flux, model) : compute_flux(sim.flux, model, sim.time.t)

    # Bundle cascading functions for future easy access
    # TODO: can we do this in a more elegant way? Maybe we should have a dedicated struct?
    # Or maybe we should simply put them in DegradationCache?
    cascading_neutrals = (cascading_N2, cascading_O2, cascading_O)

    # Initialize solution arrays
    I0 = zeros(length(z) * length(μ_center), n_E)
    Ie = zeros(length(z) * length(μ_center), n_t, n_E)
    n_t_save = isnothing(sim.time) ? 1 : length(0:sim.time.dt_requested:t_loop[end])
    Ie_save = zeros(length(z) * length(μ_center), n_t_save, n_E)

    return TransportCache(solver, degradation, matrices, Ie, Ie_save, I0, Ie_top,
                          t_loop, phase_fcn_neutrals, B2B_fragment, cascading_neutrals)
end

"""
    initialize!(sim::AuroraSimulation)

Allocate or re-allocate the working cache for `sim`.

This step performs the expensive setup that depends on the model geometry and the
resolved time grid, but does not write any output files.
"""
function initialize!(sim::AuroraSimulation)
    sim.cache = build_transport_cache(sim)
    return sim
end

function solve!(sim::AuroraSimulation)
    if isnothing(sim.time)
        solve_steady_state!(sim)
    else
        solve_time_dependent!(sim)
    end
    return sim
end

function solve_time_dependent!(sim::AuroraSimulation)
    cache = require_cache(sim)
    time = sim.time

    fill!(cache.I0, 0.0)
    fill!(cache.Ie, 0.0)
    fill!(cache.Ie_save, 0.0)

    for i_loop in 1:time.n_loop
        fill!(cache.Ie, 0.0)
        fill!(cache.matrices.Q, 0.0)

        i_start = 1 + (i_loop - 1) * (time.n_t_per_loop - 1)
        i_stop = time.n_t_per_loop + (i_loop - 1) * (time.n_t_per_loop - 1)
        Ie_top_local = @view cache.Ie_top[:, i_start:i_stop, :]

        energy_loop!(sim, Ie_top_local, i_loop, time.n_loop)

        cache.I0 .= cache.Ie[:, end, :]
        cache.Ie_save .= @view cache.Ie[:, 1:time.CFL_factor:end, :]
        save_results(sim, cache.Ie_save, cache.t_loop, cache.I0, i_loop, time.CFL_factor)
    end

    return sim
end

function solve_steady_state!(sim::AuroraSimulation)
    cache = require_cache(sim)

    fill!(cache.I0, 0.0)
    fill!(cache.Ie, 0.0)
    fill!(cache.matrices.Q, 0.0)

    Ie_top_local = @view cache.Ie_top[:, 1, :]
    energy_loop!(sim, Ie_top_local, 1, 1)
    save_results(sim, cache.Ie, cache.t_loop, cache.I0, 1, 1)

    return sim
end

function energy_loop!(sim::AuroraSimulation, Ie_top_local, i_loop::Int, n_loop::Int)
    cache = require_cache(sim)
    model = sim.model
    E_centers = model.energy_grid.E_centers
    n_E = model.energy_grid.n

    progress = isnothing(sim.time) ?
               Progress(n_E, desc="Calculating flux") :
               Progress(n_E; desc=string("Calculating flux for loop ", i_loop, "/", n_loop), dt=1.0)

    for iE in n_E:-1:1
        B2B_inelastic_neutrals = update_matrices!(cache.matrices, model,
                                                  cache.phase_fcn_neutrals, iE,
                                                  cache.B2B_fragment)

        solve_energy_step!(sim, iE, Ie_top_local; first_iteration=(iE == n_E))

        update_Q!(cache.matrices, cache.Ie, model, cache.t_loop,
                  B2B_inelastic_neutrals, cache.cascading_neutrals,
                  iE, cache.degradation)

        next!(progress)
    end

    return sim
end

function solve_energy_step!(sim::AuroraSimulation, iE::Int, Ie_top_local; first_iteration=false)
    cache = require_cache(sim)
    model = sim.model

    if isnothing(sim.time)
        @views steady_state_scheme_optimized!(cache.Ie[:, 1, iE], model,
                                              cache.matrices, iE,
                                              Ie_top_local[:, iE], cache.solver;
                                              first_iteration)
    else
        @views Crank_Nicolson_optimized!(cache.Ie[:, :, iE], cache.t_loop, model,
                                         v_of_E(model.energy_grid.E_centers[iE]),
                                         cache.matrices, iE,
                                         Ie_top_local[:, :, iE], cache.I0[:, iE],
                                         cache.solver; first_iteration)
    end

    return sim
end
