using ProgressMeter: Progress, next!

"""
    run!(sim::AuroraSimulation)

Execute the simulation described by `sim`.

Automatically calls [`initialize!`](@ref) when needed, then dispatches to the
time-dependent or steady-state execution path.

# Examples
```julia
sim = AuroraSimulation(model, flux, 0.5, 0.001, savedir)
run!(sim)
```
"""
function run!(sim::AuroraSimulation)
    if sim.cache === nothing
        @info "Simulation not initialized, calling initialize!(sim)..."
        initialize!(sim)
    end
    if sim.save_input_flux
        save_Ie_top(sim, sim.cache.Ie_top, isnothing(sim.time) ? (1:1:1) : sim.time.t)
    end
    save_parameters(sim)
    save_neutrals(sim)
    solve!(sim)
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
    cache = get_cache(sim)
    time = sim.time

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

        energy_loop!(sim, Ie_top_local, i_loop, time.n_loop)

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

function solve_steady_state!(sim::AuroraSimulation)
    cache = get_cache(sim)

    fill!(cache.I0, 0.0)
    fill!(cache.Ie, 0.0)
    fill!(cache.matrices.Q, 0.0)

    Ie_top_local = @view cache.Ie_top[:, 1, :]
    energy_loop!(sim, Ie_top_local, 1, 1)
    save_results(sim, cache.Ie, cache.t_loop, cache.I0, 1, 1)

    return sim
end

function energy_loop!(sim::AuroraSimulation, Ie_top_local, i_loop::Int, n_loop::Int)
    cache = get_cache(sim)
    model = sim.model
    E_centers = model.energy_grid.E_centers
    n_E = model.energy_grid.n

    progress = isnothing(sim.time) ?
               Progress(n_E, desc="Calculating flux",  dt=1.0) :
               Progress(n_E; desc=string("Calculating flux for loop ", i_loop, "/", n_loop), dt=1.0)

    # Energy loop: solve transport in descending energy order.
    # High-to-low ensures cascading sources are available when solving lower energies.
    for iE in n_E:-1:1
        # Update transport matrices with current energy's scattering geometry
        B2B_inelastic_neutrals = update_matrices!(cache.matrices, model,
                                                  cache.phase_fcn_neutrals, iE,
                                                  cache.B2B_fragment)

        # Solve transport equation for current energy
        solve_energy_step!(sim, iE, Ie_top_local; first_iteration=(iE == n_E))

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

function solve_energy_step!(sim::AuroraSimulation, iE::Int, Ie_top_local; first_iteration=false)
    cache = get_cache(sim)
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

function get_cache(sim::AuroraSimulation)
    cache = sim.cache
    cache === nothing && error("Simulation not initialized. Call initialize!(sim) or run!(sim).")
    return cache
end
