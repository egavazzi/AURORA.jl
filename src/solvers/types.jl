"""
    AbstractSolver

Abstract supertype for all solver strategies.

Concrete subtypes:
- [`SteadyStateSolver`](@ref): solves each time step independently as a steady-state problem
- [`TimeDependentSolver`](@ref): time-stepping with a Crank-Nicolson scheme

The solver type controls both the numerical scheme used at each energy step and
how time parameters (if any) are resolved into a time grid.
"""
abstract type AbstractSolver end

"""
    SteadyStateSolver()
    SteadyStateSolver(t_total, dt)

Solve each time step independently as a steady-state (time-independent) transport problem.

When called **without arguments**, a single steady-state solve is performed — no time grid
is created and the input flux must use [`ConstantModulation`](@ref).

When called **with `t_total` and `dt`**, a [`UniformTimeGrid`](@ref) is built and each
point is solved independently. This enables, e.g., a flickering input where each time step
is solved in steady-state.

# Examples
```julia
# Single-step steady-state
sim = AuroraSimulation(model, flux, savedir; solver=SteadyStateSolver())

# Multi-step steady-state (each time point solved independently)
sim = AuroraSimulation(model, flux, savedir; solver=SteadyStateSolver(0.5, 0.01))
```
"""
struct SteadyStateSolver <: AbstractSolver
    t_total::Union{Nothing, Float64}
    dt::Union{Nothing, Float64}

    SteadyStateSolver() = new(nothing, nothing)

    function SteadyStateSolver(t_total, dt)
        t_total > 0 || error("t_total must be positive, got $t_total")
        dt > 0 || error("dt must be positive, got $dt")
        dt <= t_total || error("dt ($dt) must be ≤ t_total ($t_total)")
        return new(Float64(t_total), Float64(dt))
    end
end

"""
    is_multi_step(solver::SteadyStateSolver)

Return `true` if the solver is configured for multi-step steady-state.
"""
is_multi_step(solver::SteadyStateSolver) = !isnothing(solver.t_total)

"""
    TimeDependentSolver(t_total, dt; CFL_number=64, max_memory_gb=8.0, n_loop=nothing)

Solve the transport equation in a time-dependent manner with a Crank-Nicolson scheme.

# Arguments
- `t_total`: total simulation time (s)
- `dt`: time step for saving data (s). The internal time step may be finer to satisfy
  the CFL condition.

# Keyword Arguments
- `CFL_number = 64`: CFL stability factor. Crank-Nicolson is unconditionally stable,
  so large values are acceptable (it will affects the accuracy though).
- `max_memory_gb = 8.0`: memory budget (GB) used to auto-split the simulation into loops.
- `n_loop = nothing`: manually specify the number of loops. Overrides `max_memory_gb`.

# Examples
```julia
sim = AuroraSimulation(model, flux, savedir;
                       solver=TimeDependentSolver(0.5, 0.001; CFL_number=128))
```
"""
struct TimeDependentSolver <: AbstractSolver
    t_total::Float64
    dt::Float64
    CFL_number::Float64
    max_memory_gb::Float64
    n_loop::Union{Nothing, Int}
end

function TimeDependentSolver(t_total, dt; CFL_number=64, max_memory_gb=8.0, n_loop=nothing)
    t_total > 0 || error("t_total must be positive, got $t_total")
    dt > 0 || error("dt must be positive, got $dt")
    dt <= t_total || error("dt ($dt) must be ≤ t_total ($t_total)")
    CFL_number > 0 || error("CFL_number must be positive, got $CFL_number")
    max_memory_gb > 0 || error("max_memory_gb must be positive, got $max_memory_gb")
    if !isnothing(n_loop)
        n_loop > 0 || error("n_loop must be positive, got $n_loop")
    end
    return TimeDependentSolver(Float64(t_total), Float64(dt), Float64(CFL_number),
                               Float64(max_memory_gb), n_loop)
end

function Base.show(io::IO, ::SteadyStateSolver)
    print(io, "SteadyStateSolver()")
end
function Base.show(io::IO, ::MIME"text/plain", s::SteadyStateSolver)
    if is_multi_step(s)
        println(io, "SteadyStateSolver (multi-step):")
        println(io, "├── t_total: ", s.t_total, " s")
        print(io,   "└── dt:      ", s.dt, " s")
    else
        print(io, "SteadyStateSolver (single-step)")
    end
end

function Base.show(io::IO, s::TimeDependentSolver)
    print(io, "TimeDependentSolver(t_total=$(s.t_total), dt=$(s.dt))")
end
function Base.show(io::IO, ::MIME"text/plain", s::TimeDependentSolver)
    println(io, "TimeDependentSolver:")
    println(io, "├── t_total:      ", s.t_total, " s")
    println(io, "├── dt:           ", s.dt, " s")
    println(io, "├── CFL number:   ", s.CFL_number)
    println(io, "├── Memory limit: ", s.max_memory_gb, " GB")
    print(io,   "└── n_loop:       ", isnothing(s.n_loop) ? "auto" : s.n_loop)
end
