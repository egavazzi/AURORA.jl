"""
    AbstractMode

Abstract supertype for simulation modes.

Concrete subtypes:
- [`SteadyStateMode`](@ref): solves each time step independently as a steady-state problem
- [`TimeDependentMode`](@ref): time-stepping with a Crank-Nicolson scheme

The mode controls both the numerical scheme used at each energy step and how
time parameters (if any) are resolved into an internal time configuration.
"""
abstract type AbstractMode end

"""
    SteadyStateMode()
    SteadyStateMode(duration, dt)

Solve each time step independently as a steady-state (time-independent) transport problem.

When called **without arguments**, a single steady-state solve is performed — no time grid
is created and the input flux must use [`ConstantModulation`](@ref).

When called **with `duration` and `dt`**, a [`UniformTimeGrid`](@ref) is built and each
point is solved independently. This enables, e.g., a flickering input where each time step
is solved in steady-state.

# Examples
```julia
# Single-step steady-state
sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode())

# Multi-step steady-state (each time point solved independently)
sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode(0.5, 0.01))
```
"""
struct SteadyStateMode <: AbstractMode
    duration::Union{Nothing, Float64}
    dt::Union{Nothing, Float64}

    SteadyStateMode() = new(nothing, nothing)

    function SteadyStateMode(duration, dt)
        duration > 0 || error("duration must be positive, got $duration")
        dt > 0 || error("dt must be positive, got $dt")
        dt <= duration || error("dt ($dt) must be ≤ duration ($duration)")
        return new(Float64(duration), Float64(dt))
    end
end

"""
    is_multi_step(mode::SteadyStateMode)

Return `true` if the mode is configured for multi-step steady-state.
"""
is_multi_step(mode::SteadyStateMode) = !isnothing(mode.duration)

"""
    SteadyState()
    SteadyState(duration, dt)

Convenience alias for [`SteadyStateMode`](@ref).
"""
const SteadyState = SteadyStateMode

"""
    TimeDependentMode(duration, dt; CFL_number=64, max_memory_gb=8.0, n_loop=nothing)

Solve the transport equation in a time-dependent manner with a Crank-Nicolson scheme.

# Arguments
- `duration`: total simulation time (s)
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
                       mode=TimeDependentMode(0.5, 0.001; CFL_number=128))
```
"""
struct TimeDependentMode <: AbstractMode
    duration::Float64
    dt::Float64
    CFL_number::Float64
    max_memory_gb::Float64
    n_loop::Union{Nothing, Int}
end

function TimeDependentMode(duration, dt; CFL_number=64, max_memory_gb=8.0, n_loop=nothing)
    duration > 0 || error("duration must be positive, got $duration")
    dt > 0 || error("dt must be positive, got $dt")
    dt <= duration || error("dt ($dt) must be ≤ duration ($duration)")
    n_intervals = Float64(duration) / Float64(dt)
    isapprox(n_intervals, round(n_intervals); atol=1e-10, rtol=1e-10) ||
        error("duration ($duration) must be an integer multiple of dt ($dt)")
    CFL_number > 0 || error("CFL_number must be positive, got $CFL_number")
    max_memory_gb > 0 || error("max_memory_gb must be positive, got $max_memory_gb")
    if !isnothing(n_loop)
        n_loop > 0 || error("n_loop must be positive, got $n_loop")
    end
    return TimeDependentMode(Float64(duration), Float64(dt), Float64(CFL_number),
                             Float64(max_memory_gb), n_loop)
end

"""
    TimeDependent(duration, dt; CFL_number=64, max_memory_gb=8.0, n_loop=nothing)

Convenience alias for [`TimeDependentMode`](@ref).
"""
const TimeDependent = TimeDependentMode

function Base.show(io::IO, s::SteadyStateMode)
    if is_multi_step(s)
        print(io, "SteadyStateMode(duration=$(s.duration), dt=$(s.dt))")
    else
        print(io, "SteadyStateMode()")
    end
end
function Base.show(io::IO, ::MIME"text/plain", s::SteadyStateMode)
    if is_multi_step(s)
        println(io, "SteadyStateMode (multi-step):")
        println(io, "├── duration: ", s.duration, " s")
        print(io,   "└── dt:      ", s.dt, " s")
    else
        print(io, "SteadyStateMode (single-step)")
    end
end

function Base.show(io::IO, s::TimeDependentMode)
    print(io, "TimeDependentMode(duration=$(s.duration), dt=$(s.dt))")
end
function Base.show(io::IO, ::MIME"text/plain", s::TimeDependentMode)
    println(io, "TimeDependentMode:")
    println(io, "├── duration:     ", s.duration, " s")
    println(io, "├── dt:           ", s.dt, " s")
    println(io, "├── CFL number:   ", s.CFL_number)
    println(io, "├── Memory limit: ", s.max_memory_gb, " GB")
    print(io,   "└── n_loop:       ", isnothing(s.n_loop) ? "auto" : s.n_loop)
end
