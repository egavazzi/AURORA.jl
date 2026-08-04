using ProgressMeter: Progress, next!
using NCDatasets: NCDataset

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
function run!(sim::AuroraSimulation; verbose::Bool = true)
    out     = sim.output
    savedir = out.savedir

    nc_path = joinpath(savedir, "simulation_data.nc")
    if isfile(nc_path) && !out.overwrite
        error("simulation_data.nc already exists in \"$savedir\". " *
              "Pass overwrite=true to AuroraOutputManager to allow overwriting.")
    end

    analysis_dir = joinpath(savedir, "analysis")
    if isfile(nc_path) && out.overwrite && isdir(analysis_dir)
        rm(analysis_dir; recursive=true)
    end

    if needs_initialization(sim)
        initialize!(sim; verbose)
    end

    mkpath(savedir)
    mkpath(joinpath(savedir, "inputs"))

    write_config_toml(sim)
    write_atmosphere_nc(sim)
    write_physics_jld2(sim)

    ds = create_simulation_nc(sim)
    try
        solve!(sim, ds; verbose)
    finally
        close(ds)
    end
    return sim
end

solve!(sim::AuroraSimulation, ds::NCDataset; verbose::Bool = true) =
    solve!(sim, sim.mode, ds; verbose)

function solve!(sim::AuroraSimulation, ::TimeDependentMode, ds::NCDataset; verbose::Bool = true)
    verbose && @info "Starting time-dependent simulation..."
    workspace = get_workspace(sim)
    time = sim.time::RefinedTimeGrid
    n_E = sim.model.energy_grid.n

    fill!(workspace.I0, 0.0)
    fill!(workspace.Ie_save, 0.0)

    for i_loop in 1:time.n_loop
        fill!(workspace.Ie, 0.0)
        fill!(workspace.matrices.Q, 0.0)

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
        Ie_top_local = @view workspace.Ie_top[:, start_internal_loop : stop_internal_loop, :]
        # Actual internal time range for this loop (may be shorter for the last loop)
        t_current_loop = workspace.t_loop[1:n_internal_loop]

        progress = Progress(n_E;
                            desc = string("Solving [loop ", i_loop, "/", time.n_loop, "]"),
                            dt = 1.0)

        energy_loop!(sim, Ie_top_local, t_current_loop, progress)

        # Save current loop final state to I0 for continuity to next loop.
        # Use n_internal_loop (not `end`) because the last loop can be shorter (Ie is
        # pre-allocated in the workspace to the maximum size n_t_per_loop).
        workspace.I0 .= @view workspace.Ie[:, n_internal_loop, :]

        # Subsample output: keep every CFL_factor-th column (gives n_save_loop + 1 columns,
        # the +1 being the shared boundary/I0 point between loops).
        n_save_cols = n_save_loop + 1
        workspace.Ie_save[:, 1:n_save_cols, :] .=
            @view workspace.Ie[:, 1:time.CFL_factor:n_internal_loop, :]

        # Build the save-time vector and Ie view for the current loop/file.
        # Loop 1: include the t=0 boundary point.
        # Loops 2+: skip the boundary (already saved in the previous file, no overlap).
        skip           = i_loop > 1     # 0 for first loop, 1 for all others
        t_save_tofile  = time.t_save[start_save_loop + skip : stop_save_loop]
        Ie_save_tofile = @view workspace.Ie_save[:, 1 + skip : n_save_cols, :]

        append_chunk_nc!(ds, Ie_save_tofile, t_save_tofile, sim)
    end

    return sim
end

solve!(sim::AuroraSimulation, ::SteadyStateMode, ds::NCDataset; verbose::Bool = true) =
    solve!(sim, sim.time, ds; verbose)
solve!(sim::AuroraSimulation, ::SingleStepConfig, ds::NCDataset; verbose::Bool = true) =
    solve_single_step_steady_state!(sim, ds; verbose)
solve!(sim::AuroraSimulation, ::UniformTimeGrid, ds::NCDataset; verbose::Bool = true) =
    solve_multi_step_steady_state!(sim, ds; verbose)

function solve_single_step_steady_state!(sim::AuroraSimulation, ds::NCDataset; verbose::Bool = true)
    verbose && @info "Starting single-step steady-state simulation..."
    workspace = get_workspace(sim)
    n_E = sim.model.energy_grid.n

    fill!(workspace.I0, 0.0)
    fill!(workspace.Ie, 0.0)
    fill!(workspace.matrices.Q, 0.0)

    Ie_top_local = @view workspace.Ie_top[:, 1, :]
    progress = Progress(n_E; desc="Solving", dt=1.0)
    energy_loop!(sim, Ie_top_local, workspace.t_loop, progress)
    append_chunk_nc!(ds, workspace.Ie[:, 1:1, :], [0.0], sim)

    return sim
end

function solve_multi_step_steady_state!(sim::AuroraSimulation, ds::NCDataset; verbose::Bool = true)
    verbose && @info "Starting multi-step steady-state simulation..."
    workspace = get_workspace(sim)
    time = sim.time::UniformTimeGrid
    n_E = sim.model.energy_grid.n
    progress = Progress(time.n_steps * n_E; desc="Solving", dt=1.0)

    for i_step in 1:time.n_steps
        # Reset state for each independent step
        fill!(workspace.I0, 0.0)
        fill!(workspace.Ie, 0.0)
        fill!(workspace.matrices.Q, 0.0)

        # Extract the boundary condition for this time step
        Ie_top_local = @view workspace.Ie_top[:, i_step, :]
        energy_loop!(sim, Ie_top_local, workspace.t_loop, progress)

        append_chunk_nc!(ds, workspace.Ie[:, 1:1, :], [time.t[i_step]], sim)
    end

    return sim
end

function energy_loop!(sim::AuroraSimulation, Ie_top_local, t_loop, progress::Progress)
    workspace = get_workspace(sim)
    model = sim.model
    n_E = model.energy_grid.n

    # Energy loop: solve transport in descending energy order.
    # High-to-low ensures cascading sources are available when solving lower energies.
    for iE in n_E:-1:1
        # Update transport matrices with current energy's scattering geometry
        B2B_inelastic_neutrals =
            update_matrices!(workspace.matrices, model, iE, workspace.B2B_kernel)

        # Solve transport equation for current energy
        solve_energy_step!(sim, sim.mode, iE, Ie_top_local, t_loop)

        # Update source term Q for lower energies from current energy's:
        # - inelastic scattering collisions → degradation → lower energies
        # - ionization collisions → cascading secondaries & degraded primaries
        update_Q!(workspace.matrices, workspace.Ie, model, workspace.t_loop,
                  B2B_inelastic_neutrals, iE, workspace.degradation)

        next!(progress)
    end

    return sim
end

function solve_energy_step!(sim::AuroraSimulation, ::SteadyStateMode,
                             iE::Int, Ie_top_local, t_loop)
    workspace = get_workspace(sim)
    model = sim.model

    @views steady_state_scheme!(workspace.Ie[:, 1, iE], model,
                                workspace.matrices, iE,
                                Ie_top_local[:, iE], workspace.solver)
    return sim
end

function solve_energy_step!(sim::AuroraSimulation, ::TimeDependentMode,
                             iE::Int, Ie_top_local, t_loop)
    workspace = get_workspace(sim)
    model = sim.model
    n_int = length(t_loop)

    @views Crank_Nicolson!(workspace.Ie[:, 1:n_int, iE], t_loop, model,
                           v_of_E(model.energy_grid.E_centers[iE]),
                           workspace.matrices, iE,
                           Ie_top_local[:, :, iE], workspace.I0[:, iE],
                           workspace.solver)
    return sim
end

function get_workspace(sim::AuroraSimulation)
    workspace = sim.workspace
    needs_initialization(sim) &&
        error("Simulation not initialized. Call initialize!(sim) or run!(sim).")
    return workspace
end
