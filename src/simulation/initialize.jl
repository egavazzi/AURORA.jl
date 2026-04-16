
"""
    initialize!(sim::AuroraSimulation)

Allocate or re-allocate the working cache for `sim`.

This step performs the expensive setup that depends on the model geometry and the
resolved time grid, but does not write any output files.
"""
function initialize!(sim::AuroraSimulation)
    @info "Initializing simulation..."
    sim.cache = build_simulation_cache(sim)
    return nothing
end

function build_simulation_cache(sim::AuroraSimulation)
    # Extract model geometry and grids
    model = sim.model
    z = model.altitude_grid.h
    μ_center = model.pitch_angle_grid.μ_center
    neutral_densities = n_neutrals(model.ionosphere)
    n_E = model.energy_grid.n

    # Set up time grid dimensions for working arrays
    n_t, t_loop = _cache_time_params(sim)

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
    Ie_top = _compute_input_flux(sim)

    # Bundle cascading functions for future easy access
    # TODO: can we do this in a more elegant way? Maybe we should have a dedicated struct?
    # Or maybe we should simply put them in DegradationCache?
    cascading_neutrals = (cascading_N2, cascading_O2, cascading_O)
    # Pre-load/calculate the cascading transfer matrices.
    preload_cascading_matrices!(model, cascading_neutrals)

    # Initialize solution arrays
    I0 = zeros(length(z) * length(μ_center), n_E)
    Ie = zeros(length(z) * length(μ_center), n_t, n_E)
    n_t_save = _save_time_length(sim, t_loop)
    Ie_save = zeros(length(z) * length(μ_center), n_t_save, n_E)

    return SimulationCache(solver, degradation, matrices, Ie, Ie_save, I0, Ie_top,
                           t_loop, phase_fcn_neutrals, B2B_fragment, cascading_neutrals)
end

# Time parameters for working arrays: (n_t, t_loop)
_cache_time_params(sim::AuroraSimulation) = _cache_time_params(sim, sim.solver)

function _cache_time_params(sim::AuroraSimulation, ::SteadyStateSolver)
    # Both single-step and multi-step SS solve one time slice at a time
    return 1, (1:1:1)
end
function _cache_time_params(sim::AuroraSimulation, ::TimeDependentSolver)
    time = sim.time::RefinedTimeGrid
    n_t = time.n_t_per_loop
    t_loop = range(0.0, step=time.dt_resolved, length=time.n_t_per_loop)
    return n_t, t_loop
end

# Compute input flux: dispatch on solver type
_compute_input_flux(sim::AuroraSimulation) = _compute_input_flux(sim, sim.solver)

function _compute_input_flux(sim::AuroraSimulation, solver::SteadyStateSolver)
    if is_multi_step(solver)
        # Multi-step SS: compute flux at all time points (like TD)
        return compute_flux(sim.flux, sim.model, sim.time.t)
    else
        # Single-step SS: no time grid
        return compute_flux(sim.flux, sim.model)
    end
end
function _compute_input_flux(sim::AuroraSimulation, ::TimeDependentSolver)
    return compute_flux(sim.flux, sim.model, sim.time.t)
end

# Number of time steps to save
_save_time_length(sim::AuroraSimulation, t_loop) = _save_time_length(sim, sim.solver, t_loop)

function _save_time_length(sim::AuroraSimulation, solver::SteadyStateSolver, t_loop)
    # Multi-step SS: save all steps; single-step SS: save 1
    return is_multi_step(solver) ? sim.time.n_steps : 1
end
function _save_time_length(sim::AuroraSimulation, ::TimeDependentSolver, t_loop)
    time = sim.time::RefinedTimeGrid
    return length(0:time.dt_requested:t_loop[end])
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

# Trigger the first (potentially expensive) load-or-compute step for all three species
# cascading functions. The closures store their transfer matrices in module-level `let`
# blocks indexed by E_grid; once loaded they stay warm for the session.
function preload_cascading_matrices!(model::AuroraModel, cascading_functions::Tuple)
    E_grid = model.energy_grid.E_edges[1:end-1]
    dE = model.energy_grid.ΔE
    # Random ionization threshold (non-physical)
    E_threshold = 15.0
    # Any primary energy above the ionization threshold works
    E_primary = 40.0
    # We just use some random numbers for the ionization thresholds
    for cascading_func in cascading_functions
        cascading_func(E_grid, dE, E_primary, E_threshold, "s")
    end
    return nothing
end
