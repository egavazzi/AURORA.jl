
"""
    initialize!(sim::AuroraSimulation;
                force_recompute=false,
                save_cache=true,
                cache_root=default_cache_root())

Allocate or re-allocate the working cache for `sim`.

# Keywords
- `force_recompute`: ignore compatible on-disk cascading caches and rebuild them
- `save_cache`: if `false`, skip writing newly computed cascading caches to disk
- `cache_root`: parent directory that contains the cascading and scattering cache folders
"""
function initialize!(sim::AuroraSimulation;
                     force_recompute::Bool = false,
                     save_cache::Bool = true,
                     cache_root::String = default_cache_root(),
                     verbose::Bool = true)
    verbose && @info "Initializing simulation..."
    cache_policy = CachePolicy(; force_recompute, save_cache, cache_root)
    if !sim.model.initialized
        initialize!(sim.model; policy=cache_policy, verbose)
    end
    # Rebuild the time configuration from the (possibly changed) model grids.
    sim.time = build_time_config(sim.model, sim.mode; verbose)
    sim.cache = build_simulation_cache(sim; cache_policy)
    sim.cache_initialized = true
    return nothing
end

needs_initialization(sim::AuroraSimulation) = !sim.cache_initialized || !sim.model.initialized

function build_dummy_simulation_cache(model::AuroraModel, time::AbstractTimeConfig)
    N_neutrals = length(model.species)
    _, t_loop = get_time_parameters(time)

    solver = SolverCache()
    degradation = DegradationCache{N_neutrals}(1, 1, 1, 1)
    matrices = TransportMatrices(1, 1, 1, 1)
    Ie = zeros(1, 1, 1)
    Ie_save = zeros(1, 1, 1)
    I0 = zeros(1, 1)
    Ie_top = zeros(1, 1, 1)
    B2B_fragment = zeros(1, 1, 1)

    return SimulationCache(solver, degradation, matrices, Ie, Ie_save, I0,
                           Ie_top, t_loop, B2B_fragment)
end

function build_simulation_cache(sim::AuroraSimulation; cache_policy::CachePolicy = CachePolicy())
    # Extract model geometry and grids
    model = sim.model
    z = model.altitude_grid.h
    μ_center = model.pitch_angle_grid.μ_center
    N_neutrals = length(model.species)
    n_E = model.energy_grid.n

    # Set up time grid dimensions for working arrays
    n_t, t_loop = get_time_parameters(sim)

    # Initialize solver and physical process caches
    solver = SolverCache()
    degradation = DegradationCache{N_neutrals}(length(μ_center), n_t, length(z), n_E)
    matrices = initialize_transport_matrices(model, t_loop)
    update_D!(matrices.D, model)
    update_Ddiffusion!(matrices.Ddiffusion, model)
    update_Mmirror!(matrices.Mmirror, model)

    # Pre-compute beam-to-beam scattering fragment
    B2B_fragment = prepare_beams2beams(model.scattering.Ω_subbeam_relative, model.scattering.P_scatter)

    # Pre-compute input flux data
    Ie_top = compute_input_flux(sim)

    # Initialize solution arrays
    I0 = zeros(length(z) * length(μ_center), n_E)
    Ie = zeros(length(z) * length(μ_center), n_t, n_E)
    n_t_save = n_steps_to_save(sim, t_loop)
    Ie_save = zeros(length(z) * length(μ_center), n_t_save, n_E)

    return SimulationCache(solver, degradation, matrices, Ie, Ie_save, I0,
                           Ie_top, t_loop, B2B_fragment)
end

# Time parameters for working arrays: (n_t, t_loop)
get_time_parameters(sim::AuroraSimulation) = get_time_parameters(sim.time)
function get_time_parameters(::SingleStepConfig)
    # Single-step steady-state solves one time slice
    return 1, (1:1:1)
end
function get_time_parameters(::UniformTimeGrid)
    # Multi-step steady-state still solves one time slice at a time
    return 1, (1:1:1)
end
function get_time_parameters(time::RefinedTimeGrid)
    n_t = time.n_t_per_loop
    t_loop = range(0.0, step=time.dt_internal, length=time.n_t_per_loop)
    return n_t, t_loop
end

# Compute input flux: dispatch on time config type
compute_input_flux(sim::AuroraSimulation) = compute_input_flux(sim, sim.time)
compute_input_flux(sim::AuroraSimulation, ::SingleStepConfig) = compute_flux(sim.flux, sim.model)
compute_input_flux(sim::AuroraSimulation, time::UniformTimeGrid) = compute_flux(sim.flux, sim.model, time.t)
compute_input_flux(sim::AuroraSimulation, time::RefinedTimeGrid) = compute_flux(sim.flux, sim.model, time.t)

# Number of time steps to save
n_steps_to_save(sim::AuroraSimulation, t_loop) = n_steps_to_save(sim.time, t_loop)
n_steps_to_save(::SingleStepConfig, t_loop) = 1
n_steps_to_save(time::UniformTimeGrid, t_loop) = time.n_steps
# n_save_per_loop + 1 to include the boundary/I0 column at the start of each loop
n_steps_to_save(time::RefinedTimeGrid, t_loop) = time.n_save_per_loop + 1
