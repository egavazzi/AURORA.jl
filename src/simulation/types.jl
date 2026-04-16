"""
    AbstractTimeGrid

Abstract supertype for time grids used by simulations.

Concrete subtypes:
- [`RefinedTimeGrid`](@ref): CFL-refined grid for time-dependent simulations
- [`UniformTimeGrid`](@ref): simple uniform grid for multi-step steady-state simulations
"""
abstract type AbstractTimeGrid end

struct RefinedTimeGrid{T<:AbstractRange{Float64}} <: AbstractTimeGrid
    t_total::Float64
    dt_requested::Float64
    dt_resolved::Float64
    CFL_factor::Int
    t::T
    n_loop::Int
    n_t_per_loop::Int
end

function RefinedTimeGrid(model::AuroraModel, solver::TimeDependentSolver)
    t_total = solver.t_total
    dt = solver.dt
    CFL_number = solver.CFL_number
    max_memory_gb = solver.max_memory_gb
    n_loop = solver.n_loop

    # Extract grid dimensions from the model
    z = model.altitude_grid.h
    n_μ = length(model.pitch_angle_grid.μ_center)
    n_E = model.energy_grid.n
    v_max = v_of_E(maximum(model.energy_grid.E_centers))

    # Apply CFL criteria to resolve the time grid and calculate CFL factor
    t_resolved, CFL_factor = CFL_criteria(t_total, dt, z, v_max, CFL_number)
    # Determine number of loops based on memory constraints, or use provided value
    n_loop_resolved = isnothing(n_loop) ? calculate_n_loop(t_resolved, length(z), n_μ, n_E;
                                                           max_memory_gb=max_memory_gb) : Int(n_loop)
    # Check if it actually fits in RAM
    check_n_loop(n_loop_resolved, length(z), n_μ, length(t_resolved), n_E)

    # Calculate timesteps per loop and actual resolved timestep (can be different from the saving dt)
    n_t_per_loop = (length(t_resolved) - 1) ÷ n_loop_resolved + 1
    dt_resolved = length(t_resolved) > 1 ? Float64(t_resolved[2] - t_resolved[1]) : Float64(dt)

    # Return fully configured time grid
    return RefinedTimeGrid(Float64(t_total), Float64(dt), dt_resolved,
                            CFL_factor, t_resolved, n_loop_resolved, n_t_per_loop)
end

function Base.show(io::IO, time::RefinedTimeGrid)
    print(io, "RefinedTimeGrid(t_total=$(time.t_total), dt=$(time.dt_requested), n_loop=$(time.n_loop))")
end

function Base.show(io::IO, ::MIME"text/plain", time::RefinedTimeGrid)
    println(io, "RefinedTimeGrid:")
    println(io, "├── Total time:    ", time.t_total, " s")
    println(io, "├── Requested dt:  ", time.dt_requested, " s")
    println(io, "├── Resolved dt:   ", time.dt_resolved, " s")
    println(io, "├── CFL factor:    ", time.CFL_factor)
    println(io, "├── n_loop:        ", time.n_loop)
    print(io,   "└── steps / loop:  ", time.n_t_per_loop)
end

"""
    UniformTimeGrid(t_total, dt)

A simple uniform time grid for multi-step steady-state simulations.

Each time point is solved independently — no CFL refinement or loop splitting is needed.
"""
struct UniformTimeGrid{T<:AbstractRange{Float64}} <: AbstractTimeGrid
    t_total::Float64
    dt::Float64
    t::T
    n_steps::Int
end

function UniformTimeGrid(t_total, dt)
    t = range(0, Float64(t_total), step=Float64(dt))
    return UniformTimeGrid(Float64(t_total), Float64(dt), t, length(t))
end

function Base.show(io::IO, time::UniformTimeGrid)
    print(io, "UniformTimeGrid(t_total=$(time.t_total), dt=$(time.dt), n_steps=$(time.n_steps))")
end

function Base.show(io::IO, ::MIME"text/plain", time::UniformTimeGrid)
    println(io, "UniformTimeGrid:")
    println(io, "├── Total time: ", time.t_total, " s")
    println(io, "├── dt:         ", time.dt, " s")
    print(io,   "└── n_steps:    ", time.n_steps)
end

"""
    AuroraSimulation{M<:AuroraModel, F<:InputFlux, S<:AbstractSolver}

A complete simulation configuration, bundling the physical model, input flux,
solver strategy, output directory, resolved time grid, and lazy simulation cache.

Use [`initialize!`](@ref) to allocate the workspace explicitly, or call [`run!`](@ref)
directly to auto-initialize and execute the simulation.

# Examples
```julia
# Single-step steady-state
sim = AuroraSimulation(model, flux, savedir; solver=SteadyStateSolver())

# Multi-step steady-state
sim = AuroraSimulation(model, flux, savedir; solver=SteadyStateSolver(0.5, 0.01))

# Time-dependent
sim = AuroraSimulation(model, flux, savedir;
                       solver=TimeDependentSolver(0.5, 0.001; CFL_number=128))
```
"""
mutable struct AuroraSimulation{M<:AuroraModel, F<:InputFlux, S<:AbstractSolver}
    const model::M
    const flux::F
    const solver::S
    const savedir::String
    const time::Union{AbstractTimeGrid, Nothing}
    const save_input_flux::Bool
    cache::Union{Nothing, SimulationCache}
end

function AuroraSimulation(model::AuroraModel, flux::InputFlux, savedir;
                          solver::AbstractSolver=SteadyStateSolver(),
                          save_input_flux=true)
    time = _build_time_grid(model, solver)
    return AuroraSimulation(model, flux, solver, String(savedir), time, save_input_flux, nothing)
end

# Build the appropriate time grid based on the solver type
_build_time_grid(model::AuroraModel, solver::SteadyStateSolver) =
    is_multi_step(solver) ? UniformTimeGrid(solver.t_total, solver.dt) : nothing
_build_time_grid(model::AuroraModel, solver::TimeDependentSolver) =
    RefinedTimeGrid(model, solver)

function Base.show(io::IO, sim::AuroraSimulation)
    mode = sim.solver isa SteadyStateSolver ? "steady-state" : "time-dependent"
    print(io, "AuroraSimulation($mode)")
end

function Base.show(io::IO, ::MIME"text/plain", sim::AuroraSimulation)
    mode = sim.solver isa SteadyStateSolver ? "Steady-state" : "Time-dependent"
    println(io, "AuroraSimulation ($mode):")
    println(io, "├── Model:       ", sim.model)
    println(io, "├── Flux:        ", sim.flux)
    println(io, "├── Solver:      ", sim.solver)
    println(io, "├── Savedir:     ", sim.savedir)
    if sim.time isa RefinedTimeGrid
        println(io, "├── t_total:     ", sim.time.t_total, " s")
        println(io, "├── dt request:  ", sim.time.dt_requested, " s")
        println(io, "├── dt resolved: ", sim.time.dt_resolved, " s")
        println(io, "├── CFL factor:  ", sim.time.CFL_factor)
        println(io, "├── n_loop:      ", sim.time.n_loop)
    elseif sim.time isa UniformTimeGrid
        println(io, "├── t_total:     ", sim.time.t_total, " s")
        println(io, "├── dt:          ", sim.time.dt, " s")
        println(io, "├── n_steps:     ", sim.time.n_steps)
    else
        println(io, "├── Time:        single-step")
    end
    println(io, "├── Cache:       ", sim.cache === nothing ? "not initialized" : "initialized")
    print(io,   "└── Save flux:   ", sim.save_input_flux)
end
