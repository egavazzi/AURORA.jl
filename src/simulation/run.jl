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
    if !sim.cache_initialized
        initialize!(sim)
    end
    if sim.save_input_flux
        t_save = save_time_axis(sim.time)
        save_Ie_top(sim, sim.cache.Ie_top, t_save)
    end
    save_parameters(sim)
    save_neutrals(sim)
    solve!(sim)
    return sim
end

save_time_axis(::SingleStepConfig) = (1:1:1)
save_time_axis(time::UniformTimeGrid) = time.t
save_time_axis(time::RefinedTimeGrid) = time.t

solve!(sim::AuroraSimulation) = solve!(sim, sim.mode)

function solve!(sim::AuroraSimulation, ::TimeDependentMode)
    @info "Starting time-dependent simulation..."
    cache = get_cache(sim)
    time = sim.time::RefinedTimeGrid
    n_E = sim.model.energy_grid.n

    fill!(cache.I0, 0.0)
    fill!(cache.Ie_save, 0.0)

    for i_loop in 1:time.n_loop
        fill!(cache.Ie, 0.0)
        fill!(cache.matrices.Q, 0.0)

        # Determine the number of internal steps and save points for this loop.
        # All loops except the last have n_save_per_loop save intervals;
        # the last loop has the remainder (i.e. ≤ n_save_per_loop).
        n_save_loop         = loop_save_count(time, i_loop)
        start_save_loop     = loop_save_start(time, i_loop)
        stop_save_loop      = start_save_loop + n_save_loop
        n_internal_loop     = loop_internal_count(time, i_loop)
        start_internal_loop = loop_internal_start(time, i_loop)
        stop_internal_loop  = start_internal_loop + n_internal_loop - 1

        # Extract the boundary condition (input flux) for this loop window.
        # Each loop starts at the shared boundary point with the previous loop.
        Ie_top_local = @view cache.Ie_top[:, start_internal_loop : stop_internal_loop, :]
        # Actual internal time range for this loop (may be shorter for the last loop)
        t_current_loop = cache.t_loop[1:n_internal_loop]

        progress = Progress(n_E;
                            desc = string("Solving [loop ", i_loop, "/", time.n_loop, "]"),
                            dt = 1.0)

        energy_loop!(sim, Ie_top_local, t_current_loop, progress)

        # Save current loop final state to I0 for continuity to next loop.
        # Use n_internal_loop (not `end`) because last loop can be shorter (Ie from cache is
        # pre-allocated to the max size n_t_per_loop).
        cache.I0 .= @view cache.Ie[:, n_internal_loop, :]

        # Subsample output: keep every CFL_factor-th column (gives n_save_loop + 1 columns,
        # the +1 being the shared boundary/I0 point between loops).
        n_save_cols = n_save_loop + 1
        cache.Ie_save[:, 1:n_save_cols, :] .= @view cache.Ie[:, 1:time.CFL_factor:n_internal_loop, :]

        # Build the save-time vector and Ie view for the current loop/file.
        # Loop 1: include the t=0 boundary point.
        # Loops 2+: skip the boundary (already saved in the previous file, no overlap).
        skip           = i_loop > 1     # 0 for first loop, 1 for all others
        t_save_tofile  = time.t_save[start_save_loop + skip : stop_save_loop]
        Ie_save_tofile = @view cache.Ie_save[:, 1 + skip : n_save_cols, :]

        save_results(sim, Ie_save_tofile, t_save_tofile, cache.I0, i_loop)
    end

    return sim
end

solve!(sim::AuroraSimulation, ::SteadyStateMode) = solve!(sim, sim.time)
solve!(sim::AuroraSimulation, ::SingleStepConfig)  = solve_single_step_steady_state!(sim)
solve!(sim::AuroraSimulation, ::UniformTimeGrid)   = solve_multi_step_steady_state!(sim)

function solve_single_step_steady_state!(sim::AuroraSimulation)
    @info "Starting single-step steady-state simulation..."
    cache = get_cache(sim)
    n_E = sim.model.energy_grid.n

    fill!(cache.I0, 0.0)
    fill!(cache.Ie, 0.0)
    fill!(cache.matrices.Q, 0.0)

    Ie_top_local = @view cache.Ie_top[:, 1, :]
    progress = Progress(n_E; desc="Solving", dt=1.0)
    energy_loop!(sim, Ie_top_local, cache.t_loop, progress)
    save_results(sim, cache.Ie, cache.t_loop, cache.I0, 1)

    return sim
end

function solve_multi_step_steady_state!(sim::AuroraSimulation)
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
        energy_loop!(sim, Ie_top_local, cache.t_loop, progress)

        # Accumulate result into the save array
        cache.Ie_save[:, i_step, :] .= @view cache.Ie[:, 1, :]
    end

    # Save all steps to a single file
    save_results(sim, cache.Ie_save, time.t, cache.I0, 1)

    return sim
end

function energy_loop!(sim::AuroraSimulation, Ie_top_local, t_loop, progress::Progress)
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
        solve_energy_step!(sim, sim.mode, iE, Ie_top_local, t_loop)

        # Update source term Q for lower energies from current energy's:
        # - inelastic scattering collisions → degradation → lower energies
        # - ionization collisions → cascading secondaries & degraded primaries
        update_Q!(cache.matrices, cache.Ie, model, cache.t_loop,
                  B2B_inelastic_neutrals, cache.cascading,
                  iE, cache.degradation)

        next!(progress)
    end

    return sim
end

function solve_energy_step!(sim::AuroraSimulation, ::SteadyStateMode,
                             iE::Int, Ie_top_local, t_loop)
    cache = get_cache(sim)
    model = sim.model

    @views steady_state_scheme!(cache.Ie[:, 1, iE], model,
                                cache.matrices, iE,
                                Ie_top_local[:, iE], cache.solver)
    return sim
end

function solve_energy_step!(sim::AuroraSimulation, ::TimeDependentMode,
                             iE::Int, Ie_top_local, t_loop)
    cache = get_cache(sim)
    model = sim.model
    n_int = length(t_loop)

    @views Crank_Nicolson!(cache.Ie[:, 1:n_int, iE], t_loop, model,
                           v_of_E(model.energy_grid.E_centers[iE]),
                           cache.matrices, iE,
                           Ie_top_local[:, :, iE], cache.I0[:, iE],
                           cache.solver)
    return sim
end

function get_cache(sim::AuroraSimulation)
    cache = sim.cache
    !sim.cache_initialized && error("Simulation not initialized. Call initialize!(sim) or run!(sim).")
    return cache
end
