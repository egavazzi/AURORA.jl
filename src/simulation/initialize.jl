
"""
    initialize!(sim::AuroraSimulation)

Allocate or re-allocate the working cache for `sim`.

This step performs the expensive setup that depends on the model geometry and the
resolved time grid, but does not write any output files.
"""
function initialize!(sim::AuroraSimulation; force_recompute::Bool = false)
    @info "Initializing simulation..."
    sim.cache = build_simulation_cache(sim; force_recompute)
    sim.cache_initialized = true
    return nothing
end

function empty_simulation_cache(model::AuroraModel, time::AbstractTimeConfig)
    neutral_densities = n_neutrals(model.ionosphere)
    _, t_loop = _cache_time_params(time)

    solver = SolverCache()
    degradation = DegradationCache(Tuple(neutral_densities), 1, 1, 1, 1)
    cascading = CascadingCache()
    matrices = TransportMatrices(1, 1, 1, 1)
    Ie = zeros(1, 1, 1)
    Ie_save = zeros(1, 1, 1)
    I0 = zeros(1, 1)
    Ie_top = zeros(1, 1, 1)
    phase_fcn_neutrals = _empty_phase_functions(model)
    B2B_fragment = zeros(size(model.scattering.Ω_subbeam_relative, 1),
                         size(model.scattering.P_scatter, 2),
                         size(model.scattering.P_scatter, 3))

    return SimulationCache(solver, degradation, cascading, matrices, Ie, Ie_save, I0,
                           Ie_top, t_loop, phase_fcn_neutrals, B2B_fragment)
end

function _empty_phase_functions(model::AuroraModel)
    n_theta = length(model.scattering.θ_scatter)
    n_energy = model.energy_grid.n
    empty_phase_pair() = (zeros(n_theta, n_energy), zeros(n_theta, n_energy))
    return (empty_phase_pair(), empty_phase_pair(), empty_phase_pair())
end

function build_simulation_cache(sim::AuroraSimulation; force_recompute::Bool = false)
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
    degradation = DegradationCache(Tuple(neutral_densities), length(μ_center), n_t, length(z), n_E)
    matrices = initialize_transport_matrices(model, t_loop)
    update_D!(matrices.D, model)
    update_Ddiffusion!(matrices.Ddiffusion, model)

    # Pre-compute scattering
    phase_fcn_neutrals = compute_phase_functions(model)
    B2B_fragment = prepare_beams2beams(model.scattering.Ω_subbeam_relative, model.scattering.P_scatter)

    # Pre-compute input flux data
    Ie_top = _compute_input_flux(sim)

    # Build the cascading cache
    cascading = CascadingCache()
    # Pre-load/calculate the cascading transfer matrices.
    preload_cascading_matrices!(model, cascading; force_recompute)

    # Initialize solution arrays
    I0 = zeros(length(z) * length(μ_center), n_E)
    Ie = zeros(length(z) * length(μ_center), n_t, n_E)
    n_t_save = _save_time_length(sim, t_loop)
    Ie_save = zeros(length(z) * length(μ_center), n_t_save, n_E)

    return SimulationCache(solver, degradation, cascading, matrices, Ie, Ie_save, I0,
                           Ie_top, t_loop, phase_fcn_neutrals, B2B_fragment)
end

# Time parameters for working arrays: (n_t, t_loop)
_cache_time_params(sim::AuroraSimulation) = _cache_time_params(sim.time)

function _cache_time_params(::SingleStepConfig)
    # Single-step steady-state solves one time slice
    return 1, (1:1:1)
end
function _cache_time_params(::UniformTimeGrid)
    # Multi-step steady-state still solves one time slice at a time
    return 1, (1:1:1)
end
function _cache_time_params(time::RefinedTimeGrid)
    n_t = time.n_t_per_loop
    t_loop = range(0.0, step=time.dt_resolved, length=time.n_t_per_loop)
    return n_t, t_loop
end

# Compute input flux: dispatch on time config type
_compute_input_flux(sim::AuroraSimulation) = _compute_input_flux(sim, sim.time)

_compute_input_flux(sim::AuroraSimulation, ::SingleStepConfig) = compute_flux(sim.flux, sim.model)
_compute_input_flux(sim::AuroraSimulation, time::UniformTimeGrid) = compute_flux(sim.flux, sim.model, time.t)
_compute_input_flux(sim::AuroraSimulation, time::RefinedTimeGrid) = compute_flux(sim.flux, sim.model, time.t)

# Number of time steps to save
_save_time_length(sim::AuroraSimulation, t_loop) = _save_time_length(sim.time, t_loop)

_save_time_length(::SingleStepConfig, t_loop) = 1
_save_time_length(time::UniformTimeGrid, t_loop) = time.n_steps
_save_time_length(time::RefinedTimeGrid, t_loop) = length(0:time.dt_requested:t_loop[end])

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

# Trigger the first (potentially expensive) load-or-compute step for all three species.
function preload_cascading_matrices!(model::AuroraModel, cascading::CascadingCache;
                                     force_recompute::Bool = false)
    for species_cache in cascading
        ensure_cascading_loaded!(species_cache, model.energy_grid; force_recompute)
    end
    return nothing
end
