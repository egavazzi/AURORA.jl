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
        t_save = input_time_axis(sim.time)
        save_Ie_top(sim, sim.cache.Ie_top, t_save)
    end
    save_parameters(sim)
    save_neutrals(sim)
    solve!(sim)
    return sim
end

input_time_axis(::SingleStepConfig) = (1:1:1)
input_time_axis(time::UniformTimeGrid) = time.t
input_time_axis(time::RefinedTimeGrid) = time.t

save_time_axis(::SingleStepConfig) = (1:1:1)
save_time_axis(time::UniformTimeGrid) = time.t
save_time_axis(time::RefinedTimeGrid) = time.save_t

solve!(sim::AuroraSimulation) = solve!(sim, sim.mode)

function solve!(sim::AuroraSimulation, ::TimeDependentMode)
    @info "Starting time-dependent simulation..."
    cache = get_cache(sim)
    time = sim.time::RefinedTimeGrid
    n_E = sim.model.energy_grid.n

    fill!(cache.I0, 0.0)
    fill!(cache.Ie, 0.0)
    fill!(cache.Ie_save, 0.0)

    i_savefile = 0
    for (i_loop, loop_slice) in enumerate(time.loop_slices)
        fill!(cache.Ie, 0.0)
        fill!(cache.matrices.Q, 0.0)

        n_t_loop = length(loop_slice.resolved_indices)
        t_loop_local = cache.t_loop[1:n_t_loop]
        Ie_local = @view cache.Ie[:, 1:n_t_loop, :]
        Q_local = @view cache.matrices.Q[:, 1:n_t_loop, :]
        Ie_top_local = @view cache.Ie_top[:, loop_slice.resolved_indices, :]
        progress = Progress(n_E;
                    desc=string("Solving [loop ", i_loop, "/", time.n_loop, "]"),
                    dt=1.0)

        energy_loop!(sim, Ie_local, t_loop_local, Q_local, Ie_top_local, progress)

        # Save current loop final state to I0 for continuity to next loop
        cache.I0 .= Ie_local[:, end, :]

        n_save = length(loop_slice.output_indices)
        n_save == 0 && continue

        cache.Ie_save[:, 1:n_save, :] .= @view Ie_local[:, loop_slice.local_save_indices, :]
        t_save = @view save_time_axis(time)[loop_slice.output_indices]
        length(t_save) == n_save || error("save grid and saved data length diverged for loop $i_loop")

        i_savefile += 1
        save_results(sim, view(cache.Ie_save, :, 1:n_save, :), t_save, cache.I0, i_savefile)
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
    energy_loop!(sim, cache.Ie, cache.t_loop, cache.matrices.Q, Ie_top_local, progress)
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
        energy_loop!(sim, cache.Ie, cache.t_loop, cache.matrices.Q, Ie_top_local, progress)

        # Accumulate result into the save array
        cache.Ie_save[:, i_step, :] .= @view cache.Ie[:, 1, :]
    end

    # Save all steps to a single file
    save_results(sim, cache.Ie_save, time.t, cache.I0, 1)

    return sim
end

function energy_loop!(sim::AuroraSimulation, Ie, t_loop, Q, Ie_top_local, progress::Progress)
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
        solve_energy_step!(sim, sim.mode, Ie, t_loop, Q, iE, Ie_top_local)

        # Update source term Q for lower energies from current energy's:
        # - inelastic scattering collisions → degradation → lower energies
        # - ionization collisions → cascading secondaries & degraded primaries
        update_Q!(Q, Ie, model, t_loop,
                  B2B_inelastic_neutrals, cache.cascading,
                  iE, cache.degradation)

        next!(progress)
    end

    return sim
end

function solve_energy_step!(sim::AuroraSimulation, ::SteadyStateMode,
                             Ie, t_loop, Q, iE::Int, Ie_top_local)
    cache = get_cache(sim)
    model = sim.model

    @views steady_state_scheme!(Ie[:, 1, iE], model,
                                cache.matrices, Q, iE,
                                Ie_top_local[:, iE], cache.solver)
    return sim
end

function solve_energy_step!(sim::AuroraSimulation, ::TimeDependentMode,
                             Ie, t_loop, Q, iE::Int, Ie_top_local)
    cache = get_cache(sim)
    model = sim.model

    @views Crank_Nicolson!(Ie[:, :, iE], t_loop, model,
                           v_of_E(model.energy_grid.E_centers[iE]),
                           cache.matrices, Q, iE,
                           Ie_top_local[:, :, iE], cache.I0[:, iE],
                           cache.solver)
    return sim
end

function get_cache(sim::AuroraSimulation)
    cache = sim.cache
    !sim.cache_initialized && error("Simulation not initialized. Call initialize!(sim) or run!(sim).")
    return cache
end
