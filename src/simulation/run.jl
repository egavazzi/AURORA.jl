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
        check_bottom_boundary(sim, ds)
    finally
        close(ds)
    end
    return sim
end

solve!(sim::AuroraSimulation, ds::NCDataset; verbose::Bool = true) =
    solve!(sim, sim.mode, ds; verbose)

function solve!(sim::AuroraSimulation, ::TimeDependentMode, ds::NCDataset; verbose::Bool = true)
    verbose && @info "Starting time-dependent simulation..."
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
    cache = get_cache(sim)
    n_E = sim.model.energy_grid.n

    fill!(cache.I0, 0.0)
    fill!(cache.Ie, 0.0)
    fill!(cache.matrices.Q, 0.0)

    Ie_top_local = @view cache.Ie_top[:, 1, :]
    progress = Progress(n_E; desc="Solving", dt=1.0)
    energy_loop!(sim, Ie_top_local, cache.t_loop, progress)
    append_chunk_nc!(ds, cache.Ie[:, 1:1, :], [0.0], sim)

    return sim
end

function solve_multi_step_steady_state!(sim::AuroraSimulation, ds::NCDataset; verbose::Bool = true)
    verbose && @info "Starting multi-step steady-state simulation..."
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

    append_chunk_nc!(ds, cache.Ie_save, collect(Float64, time.t), sim)

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
        B2B_inelastic_neutrals = update_matrices!(cache.matrices, model, iE, cache.B2B_fragment)

        # Solve transport equation for current energy
        solve_energy_step!(sim, sim.mode, iE, Ie_top_local, t_loop)

        # Update source term Q for lower energies from current energy's:
        # - inelastic scattering collisions → degradation → lower energies
        # - ionization collisions → cascading secondaries & degraded primaries
        update_Q!(cache.matrices, cache.Ie, model, cache.t_loop,
                  B2B_inelastic_neutrals, iE, cache.degradation)

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

"""
    check_bottom_boundary(sim::AuroraSimulation, ds; threshold=0.01)

Post-solve sanity check of the bottom boundary placement, called automatically by
[`run!`](@ref). Compares the time-averaged downward electron energy flux reaching the
bottom of the altitude grid with the one injected at the top. The bottom boundary
absorbs any remaining flux, so a non-negligible ratio means the electrons are not fully
stopped within the grid and the low-altitude part of the solution is clipped. Emits a
warning with a suggested bottom altitude (from [`suggest_bottom_altitude`](@ref)) when
the ratio exceeds `threshold`.

Reads the bottom z-slice of `Ie` from `ds` in blocks of time steps, so the cost stays
negligible compared to the solve.
"""
function check_bottom_boundary(sim::AuroraSimulation, ds::NCDataset; threshold = 0.01)
    E = sim.model.energy_grid.E_centers
    # Downward beams (μ < 0) come first in beam order (θ_lims goes 180° → 0°)
    n_down = count(<(0), sim.model.pitch_angle_grid.μ_center)

    # Time-averaged downward energy flux: sum over down beams and energy bins, mean over
    # time, read in blocks of time steps to keep memory bounded
    function mean_down_eflux(read_block, n_t)
        s = 0.0
        # Load in arbitrary blocks of 128 time steps. Reasonably small to keep memory low,
        # but large enough to amortize the overhead of reading from disk.
        for t1 in 1:128:n_t
            t2 = min(t1 + 127, n_t)
            s += sum(read_block(t1, t2) .* reshape(E, 1, 1, :))
        end
        return s / n_t
    end

    Ietop_v = ds["Ie_input"].var    # [n_μ, n_t_top, n_E]
    eflux_top = mean_down_eflux((t1, t2) -> Ietop_v[1:n_down, t1:t2, :], size(Ietop_v, 2))
    eflux_top <= 0 && return nothing

    # The bottom grid point itself is forced to 0 by the boundary condition, so the flux
    # about to be absorbed is the downward flux at the second grid point
    Ie_v = ds["Ie"].var     # [n_z, n_μ, n_t, n_E]
    eflux_bottom = mean_down_eflux((t1, t2) -> Ie_v[2, 1:n_down, t1:t2, :], size(Ie_v, 3))

    ratio = eflux_bottom / eflux_top
    if ratio > threshold
        bottom = sim.model.altitude_grid.bottom
        suggestion = suggest_bottom_altitude(sim.model.energy_grid.E_max,
                                             sim.model.ionosphere.msis_file)
        @warn "$(round(100 * ratio, sigdigits=2))% of the precipitating energy flux " *
              "reaches the bottom of the altitude grid ($(bottom) km), where it is " *
              "absorbed by the boundary condition. The solution is likely clipped at " *
              "low altitudes. Consider extending the grid down to ~$(suggestion) km " *
              "(see suggest_bottom_altitude)."
    end
    return nothing
end
