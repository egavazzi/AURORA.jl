using ProgressMeter: Progress, next!

"""
    run!(sim::AuroraSimulation)

Execute the simulation described by `sim`.

Automatically calls [`initialize!`](@ref) when needed, then dispatches to the
appropriate execution path based on the selected mode.

# Examples
```julia
sim = AuroraSimulation(model, flux, savedir;
                       mode=TimeDependentMode(duration=0.5, dt=0.001, CFL_number=128))
run!(sim)
```
"""
function run!(sim::AuroraSimulation)
    if sim.cache === nothing
        initialize!(sim)
    end
    if sim.save_input_flux
        t_save = _save_time_axis(sim.time)
        save_Ie_top(sim, sim.cache.Ie_top, t_save)
    end
    save_parameters(sim)
    save_neutrals(sim)
    solve!(sim)
    return sim
end

_save_time_axis(::SingleStepConfig) = (1:1:1)
_save_time_axis(time::UniformTimeGrid) = time.t
_save_time_axis(time::RefinedTimeGrid) = time.t

solve!(sim::AuroraSimulation) = _solve!(sim, sim.mode)

function _solve!(sim::AuroraSimulation, ::TimeDependentMode)
    @info "Starting time-dependent simulation..."
    cache = get_cache(sim)
    time = sim.time::RefinedTimeGrid
    n_E = sim.model.energy_grid.n

    fill!(cache.I0, 0.0)
    fill!(cache.Ie, 0.0)
    fill!(cache.Ie_save, 0.0)

    for i_loop in 1:time.n_loop
        fill!(cache.Ie, 0.0)
        fill!(cache.matrices.Q, 0.0)

        # Extract the boundary condition (input flux) for this loop window.
        # Windows overlap by 1 timestep to ensure continuity between loops:
        # e.g. if n_t_per_loop=10, loop 1 gets timesteps [1-10], loop 2 gets [10-19],
        # loop 3 gets [19-28], etc.
        i_start = 1 + (i_loop - 1) * (time.n_t_per_loop - 1)
        i_stop = time.n_t_per_loop + (i_loop - 1) * (time.n_t_per_loop - 1)
        Ie_top_local = @view cache.Ie_top[:, i_start:i_stop, :]
        progress = Progress(n_E;
                    desc=string("Solving [loop ", i_loop, "/", time.n_loop, "]"),
                    dt=1.0)

        energy_loop!(sim, Ie_top_local, progress)

        # Save current loop final state to I0 for continuity to next loop
        cache.I0 .= cache.Ie[:, end, :]
        # Subsample output: save every CFL_factor-th timestep to disk (we don't want to
        # save for _every_ timestep resolved internally)
        cache.Ie_save .= @view cache.Ie[:, 1:time.CFL_factor:end, :]
        # And save current loop to disk
        save_results(sim, cache.Ie_save, cache.t_loop, cache.I0, i_loop, time.CFL_factor)
    end

    return sim
end

_solve!(sim::AuroraSimulation, ::SteadyStateMode) = _solve!(sim, sim.time)
_solve!(sim::AuroraSimulation, ::SingleStepConfig)  = _solve_single_step_steady_state!(sim)
_solve!(sim::AuroraSimulation, ::UniformTimeGrid)   = _solve_multi_step_steady_state!(sim)

function _solve_single_step_steady_state!(sim::AuroraSimulation)
    @info "Starting single-step steady-state simulation..."
    cache = get_cache(sim)
    n_E = sim.model.energy_grid.n

    fill!(cache.I0, 0.0)
    fill!(cache.Ie, 0.0)
    fill!(cache.matrices.Q, 0.0)

    Ie_top_local = @view cache.Ie_top[:, 1, :]
    progress = Progress(n_E; desc="Solving", dt=1.0)
    energy_loop!(sim, Ie_top_local, progress)
    save_results(sim, cache.Ie, cache.t_loop, cache.I0, 1, 1)

    return sim
end

function _solve_multi_step_steady_state!(sim::AuroraSimulation)
    @info "Starting multi-step steady-state simulation..."
    cache = get_cache(sim)
    time = sim.time::UniformTimeGrid
    n_E = sim.model.energy_grid.n
    progress = Progress(time.n_steps * n_E; desc="Solving", dt=1.0)

    for i_step in 1:time.n_steps
        # Reset state for each independent step
        fill!(cache.I0, 0.0)
        fill!(cache.Ie, 0.0)
        fill!(cache.matrices.Q, 0.0)

        # Extract the boundary condition for this time step
        Ie_top_local = @view cache.Ie_top[:, i_step, :]
        energy_loop!(sim, Ie_top_local, progress)

        # Accumulate result into the save array
        cache.Ie_save[:, i_step, :] .= @view cache.Ie[:, 1, :]
    end

    # Save all steps to a single file
    save_results(sim, cache.Ie_save, time.t, cache.I0, 1, 1)

    return sim
end

function energy_loop!(sim::AuroraSimulation, Ie_top_local, progress::Progress)
    cache = get_cache(sim)
    model = sim.model
    n_E = model.energy_grid.n

    # Energy loop: solve transport in descending energy order.
    # High-to-low ensures cascading sources are available when solving lower energies.
    for iE in n_E:-1:1
        # Update transport matrices with current energy's scattering geometry
        B2B_inelastic_neutrals = update_matrices!(cache.matrices, model,
                                                  cache.phase_fcn_neutrals, iE,
                                                  cache.B2B_fragment)

        # Solve transport equation for current energy
        _solve_energy_step!(sim, sim.mode, iE, Ie_top_local; first_iteration=(iE == n_E))

        # Update source term Q for lower energies from current energy's:
        # - inelastic scattering collisions → degradation → lower energies
        # - ionization collisions → cascading secondaries & degraded primaries
        update_Q!(cache.matrices, cache.Ie, model, cache.t_loop,
                  B2B_inelastic_neutrals, cache.cascading_neutrals,
                  iE, cache.degradation)

        next!(progress)
    end

    return sim
end

function _solve_energy_step!(sim::AuroraSimulation, ::SteadyStateMode,
                             iE::Int, Ie_top_local; first_iteration=false)
    cache = get_cache(sim)
    model = sim.model

    @views steady_state_scheme_optimized!(cache.Ie[:, 1, iE], model,
                                          cache.matrices, iE,
                                          Ie_top_local[:, iE], cache.solver;
                                          first_iteration)
    return sim
end

function _solve_energy_step!(sim::AuroraSimulation, ::TimeDependentMode,
                             iE::Int, Ie_top_local; first_iteration=false)
    cache = get_cache(sim)
    model = sim.model

    @views Crank_Nicolson_optimized!(cache.Ie[:, :, iE], cache.t_loop, model,
                                     v_of_E(model.energy_grid.E_centers[iE]),
                                     cache.matrices, iE,
                                     Ie_top_local[:, :, iE], cache.I0[:, iE],
                                     cache.solver; first_iteration)
    return sim
end

function get_cache(sim::AuroraSimulation)
    cache = sim.cache
    cache === nothing && error("Simulation not initialized. Call initialize!(sim) or run!(sim).")
    return cache
end
