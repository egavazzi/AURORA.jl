# ──────────────────────────────────────────────────────────────────────────────
# Simulation Modes
# ──────────────────────────────────────────────────────────────────────────────

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
    SteadyStateMode(; duration, dt)

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
sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode(duration=0.5, dt=0.01))
```
"""
struct SteadyStateMode <: AbstractMode
    duration::Union{Nothing, Float64}
    dt::Union{Nothing, Float64}
end

function SteadyStateMode(; duration=nothing, dt=nothing)
    if isnothing(duration) && isnothing(dt)
        return SteadyStateMode(nothing, nothing)
    end
    isnothing(duration) && error("duration is missing (dt=$dt was provided)")
    isnothing(dt) && error("dt is missing (duration=$duration was provided)")
    duration > 0 || error("duration must be positive, got $duration")
    dt > 0 || error("dt must be positive, got $dt")
    dt <= duration || error("dt ($dt) must be ≤ duration ($duration)")
    isapprox(round(duration / dt) * dt, duration; rtol=1e-10) ||
        error("duration ($duration) must be an integer multiple of dt ($dt)")
    return SteadyStateMode(Float64(duration), Float64(dt))
end

"""
    is_multi_step(mode::SteadyStateMode)

Return `true` if the mode is configured for multi-step steady-state.
"""
is_multi_step(mode::SteadyStateMode) = !isnothing(mode.duration)

"""
    SteadyState()
    SteadyState(; duration, dt)

Convenience alias for [`SteadyStateMode`](@ref).
"""
const SteadyState = SteadyStateMode

"""
    TimeDependentMode(; duration, dt, CFL_number=64, max_memory_gb=8.0, n_loop=nothing)

Solve the transport equation in a time-dependent manner with a Crank-Nicolson scheme.

# Keyword Arguments
- `duration`: total simulation time (s)
- `dt`: time step for saving data (s). The internal time step may be finer to satisfy
  the CFL condition.
- `CFL_number = 64`: CFL stability factor. Crank-Nicolson is unconditionally stable,
  so large values are acceptable (it will affect the accuracy though).
- `max_memory_gb = 8.0`: memory budget (GB) used to auto-split the simulation into loops.
- `n_loop = nothing`: manually specify the number of loops. Overrides `max_memory_gb`.

# Examples
```julia
sim = AuroraSimulation(model, flux, savedir;
                       mode=TimeDependentMode(duration=0.5, dt=0.001, CFL_number=128))
```
"""
struct TimeDependentMode <: AbstractMode
    duration::Float64
    dt::Float64
    CFL_number::Float64
    max_memory_gb::Float64
    n_loop::Union{Nothing, Int}

    function TimeDependentMode(; duration, dt, CFL_number=64, max_memory_gb=8.0, n_loop=nothing)
        # A bunch of input validation to fail early if the user config is invalid
        duration > 0 || error("duration must be positive, got $duration")
        dt > 0 || error("dt must be positive, got $dt")
        dt <= duration || error("dt ($dt) must be ≤ duration ($duration)")
        isapprox(round(duration / dt) * dt, duration; rtol=1e-10) ||
            error("duration ($duration) must be an integer multiple of dt ($dt)")
        CFL_number > 0 || error("CFL_number must be positive, got $CFL_number")
        max_memory_gb > 0 || error("max_memory_gb must be positive, got $max_memory_gb")
        if !isnothing(n_loop)
            n_loop > 0 || error("n_loop must be positive, got $n_loop")
        end
        return new(Float64(duration), Float64(dt), Float64(CFL_number),
                   Float64(max_memory_gb), n_loop)
    end
end

"""
    TimeDependent(; duration, dt, CFL_number=64, max_memory_gb=8.0, n_loop=nothing)

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













# ──────────────────────────────────────────────────────────────────────────────
# Simulation Time Configurations
# ──────────────────────────────────────────────────────────────────────────────

"""
    AbstractTimeConfig

Abstract supertype for resolved time configurations used by simulations.

Concrete subtypes:
- [`SingleStepConfig`](@ref): explicit single-step configuration for steady-state simulations
- [`UniformTimeGrid`](@ref): simple uniform grid for multi-step steady-state simulations
- [`RefinedTimeGrid`](@ref): CFL-refined grid for time-dependent simulations
"""
abstract type AbstractTimeConfig end

"""
    SingleStepConfig()

Explicit single-step time configuration.

Used for steady-state simulations without temporal evolution.
"""
struct SingleStepConfig <: AbstractTimeConfig
    t::UnitRange{Int}
    n_steps::Int
end

SingleStepConfig() = SingleStepConfig(1:1, 1)

function Base.show(io::IO, ::SingleStepConfig)
    print(io, "SingleStepConfig()")
end

function Base.show(io::IO, ::MIME"text/plain", ::SingleStepConfig)
    print(io, "SingleStepConfig (single-step)")
end

"""
    UniformTimeGrid(duration, dt)

A simple uniform time grid for multi-step steady-state simulations.

Each time point is solved independently — no CFL refinement or loop splitting is needed.
"""
struct UniformTimeGrid{T<:AbstractRange{Float64}} <: AbstractTimeConfig
    duration::Float64
    dt::Float64
    t::T
    n_steps::Int
end

function UniformTimeGrid(duration, dt)
    duration > 0 || error("duration must be positive, got $duration")
    dt > 0 || error("dt must be positive, got $dt")
    dt <= duration || error("dt ($dt) must be ≤ duration ($duration)")

    n_intervals_round = round(Int, Float64(duration) / Float64(dt))
    n_steps = n_intervals_round + 1
    t = range(0.0, stop=Float64(duration), length=n_steps)
    return UniformTimeGrid(Float64(duration), Float64(dt), t, n_steps)
end

function Base.show(io::IO, time::UniformTimeGrid)
    print(io, "UniformTimeGrid(duration=$(time.duration), dt=$(time.dt), n_steps=$(time.n_steps))")
end

function Base.show(io::IO, ::MIME"text/plain", time::UniformTimeGrid)
    println(io, "UniformTimeGrid:")
    println(io, "├── Duration:   ", time.duration, " s")
    println(io, "├── dt:         ", time.dt, " s")
    print(io,   "└── n_steps:    ", time.n_steps)
end

struct TimeLoopSlice
    resolved_indices::UnitRange{Int}
    local_save_indices::StepRange{Int, Int}
    output_indices::UnitRange{Int}
end

function build_time_loop_slices(n_resolved::Int, n_loop::Int, save_indices::StepRange{Int, Int})
    n_intervals = n_resolved - 1
    n_intervals > 0 || error("time-dependent simulations require at least two resolved time points")
    n_loop > 0 || error("n_loop must be positive, got $n_loop")
    n_loop <= n_intervals || error("n_loop ($n_loop) must be ≤ the number of resolved time intervals ($n_intervals)")

    base_intervals = fld(n_intervals, n_loop)
    extra_intervals = mod(n_intervals, n_loop)
    save_stride = step(save_indices)
    slices = Vector{TimeLoopSlice}(undef, n_loop)

    resolved_start = 1
    for i_loop in 1:n_loop
        n_local_intervals = base_intervals + (i_loop <= extra_intervals ? 1 : 0)
        resolved_stop = resolved_start + n_local_intervals

        first_owned_resolved = i_loop == 1 ? resolved_start : resolved_start + 1
        output_start = searchsortedfirst(save_indices, first_owned_resolved)
        output_stop = searchsortedlast(save_indices, resolved_stop)
        output_indices = output_start:output_stop

        local_save_indices = if isempty(output_indices)
            1:save_stride:0
        else
            first_local_save = save_indices[first(output_indices)] - resolved_start + 1
            last_local_save = save_indices[last(output_indices)] - resolved_start + 1
            first_local_save:save_stride:last_local_save
        end

        slices[i_loop] = TimeLoopSlice(resolved_start:resolved_stop,
                                       local_save_indices,
                                       output_indices)
        resolved_start = resolved_stop
    end

    validate_time_loop_slices(slices, n_resolved, length(save_indices))

    return slices
end

function validate_time_loop_slices(slices::Vector{TimeLoopSlice}, n_resolved::Int, n_save::Int)
    isempty(slices) && error("at least one loop slice is required")
    first(first(slices).resolved_indices) == 1 || error("loop partition must start at the first resolved time index")
    last(last(slices).resolved_indices) == n_resolved || error("loop partition must end at the last resolved time index")

    covered_intervals = 0
    next_output_index = 1

    for i_slice in eachindex(slices)
        slice = slices[i_slice]
        covered_intervals += length(slice.resolved_indices) - 1

        if i_slice > 1
            prev_slice = slices[i_slice - 1]
            first(slice.resolved_indices) == last(prev_slice.resolved_indices) ||
                error("adjacent loop slices must overlap by exactly one resolved time index")
        end

        length(slice.local_save_indices) == length(slice.output_indices) ||
            error("loop slice save indices are inconsistent")

        if !isempty(slice.local_save_indices)
            first(slice.local_save_indices) >= 1 || error("local save indices must start within the loop window")
            last(slice.local_save_indices) <= length(slice.resolved_indices) ||
                error("local save indices must stay within the loop window")
            first(slice.output_indices) == next_output_index ||
                error("output save-grid ownership must be contiguous across loop slices")
            next_output_index = last(slice.output_indices) + 1
        end
    end

    covered_intervals == (n_resolved - 1) || error("loop partition must cover every resolved interval exactly once")
    next_output_index == (n_save + 1) || error("loop partition must cover every saved output step exactly once")

    return nothing
end

struct RefinedTimeGrid{T<:AbstractRange{Float64}} <: AbstractTimeConfig
    duration::Float64
    dt_requested::Float64
    dt_resolved::Float64
    CFL_factor::Int
    t::T
    n_loop::Int
    save_indices::StepRange{Int, Int}
    save_t::Vector{Float64}
    loop_slices::Vector{TimeLoopSlice}
    max_loop_length::Int
    max_save_steps::Int
end

function RefinedTimeGrid(model::AuroraModel, mode::TimeDependentMode)
    duration = mode.duration
    dt = mode.dt
    CFL_number = mode.CFL_number
    max_memory_gb = mode.max_memory_gb
    n_loop = mode.n_loop

    # Extract grid dimensions from the model
    z = model.altitude_grid.h
    n_μ = length(model.pitch_angle_grid.μ_center)
    n_E = model.energy_grid.n
    v_max = v_of_E(maximum(model.energy_grid.E_centers))

    # Apply CFL criteria to resolve the time grid and calculate CFL factor
    t_resolved, CFL_factor = CFL_criteria(duration, dt, z, v_max, CFL_number)
    # Determine number of loops based on memory constraints, or use provided value
    n_loop_resolved = isnothing(n_loop) ? calculate_n_loop(t_resolved, length(z), n_μ, n_E;
                                                           max_memory_gb=max_memory_gb) : Int(n_loop)
    # Check if it actually fits in RAM
    check_n_loop(n_loop_resolved, length(z), n_μ, length(t_resolved), n_E)

    save_indices = 1:CFL_factor:length(t_resolved)
    save_t = collect(t_resolved[save_indices])
    loop_slices = build_time_loop_slices(length(t_resolved), n_loop_resolved, save_indices)
    max_loop_length = maximum(length(slice.resolved_indices) for slice in loop_slices)
    max_save_steps = maximum(length(slice.output_indices) for slice in loop_slices)

    # Calculate actual resolved timestep (can be different from the saving dt)
    dt_resolved = length(t_resolved) > 1 ? Float64(t_resolved[2] - t_resolved[1]) : Float64(dt)

    # Return fully configured time grid
    return RefinedTimeGrid(Float64(duration), Float64(dt), dt_resolved,
                            CFL_factor, t_resolved, n_loop_resolved,
                            save_indices, save_t, loop_slices,
                            max_loop_length, max_save_steps)
end

function Base.show(io::IO, time::RefinedTimeGrid)
    print(io, "RefinedTimeGrid(duration=$(time.duration), dt=$(time.dt_requested), n_loop=$(time.n_loop))")
end

function Base.show(io::IO, ::MIME"text/plain", time::RefinedTimeGrid)
    println(io, "RefinedTimeGrid:")
    println(io, "├── Duration:      ", time.duration, " s")
    println(io, "├── Requested dt:  ", time.dt_requested, " s")
    println(io, "├── Resolved dt:   ", time.dt_resolved, " s")
    println(io, "├── CFL factor:    ", time.CFL_factor)
    println(io, "├── n_loop:        ", time.n_loop)
    print(io,   "└── max steps / loop:  ", time.max_loop_length)
end









# ──────────────────────────────────────────────────────────────────────────────
# Simulation Object
# ──────────────────────────────────────────────────────────────────────────────
"""
    AuroraSimulation{M<:AuroraModel, F<:InputFlux, S<:AbstractMode}

A complete simulation configuration, bundling the physical model, input flux,
mode strategy, output directory, resolved time configuration, and lazy simulation cache.

Use [`initialize!`](@ref) to allocate the workspace explicitly, or call [`run!`](@ref)
directly to auto-initialize and execute the simulation.

# Examples
```julia
# Single-step steady-state
sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode())

# Multi-step steady-state
sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode(duration=0.5, dt=0.01))

# Time-dependent
sim = AuroraSimulation(model, flux, savedir;
                       mode=TimeDependentMode(duration=0.5, dt=0.001, CFL_number=128))
```
"""
mutable struct AuroraSimulation{M<:AuroraModel, F<:InputFlux, S<:AbstractMode,
                               T<:AbstractTimeConfig, C<:SimulationCache}
    const model::M
    const flux::F
    const mode::S
    const savedir::String
    const time::T
    const save_input_flux::Bool
    cache::C
    cache_initialized::Bool
end

function AuroraSimulation(model::AuroraModel, flux::InputFlux, savedir;
                          mode::AbstractMode=SteadyStateMode(),
                          save_input_flux=true)
    time = _build_time_config(model, mode)
    cache = build_dummy_simulation_cache(model, time)
    return AuroraSimulation(model, flux, mode, String(savedir), time, save_input_flux,
                            cache, false)
end

# Build the appropriate time configuration based on the mode
_build_time_config(model::AuroraModel, mode::SteadyStateMode) =
    is_multi_step(mode) ? UniformTimeGrid(mode.duration, mode.dt) : SingleStepConfig()
_build_time_config(model::AuroraModel, mode::TimeDependentMode) =
    RefinedTimeGrid(model, mode)

function Base.show(io::IO, sim::AuroraSimulation)
    mode_name = sim.mode isa SteadyStateMode ? "steady-state" : "time-dependent"
    print(io, "AuroraSimulation($mode_name)")
end

function Base.show(io::IO, ::MIME"text/plain", sim::AuroraSimulation)
    mode_name = sim.mode isa SteadyStateMode ? "Steady-state" : "Time-dependent"
    println(io, "AuroraSimulation ($mode_name):")
    println(io, "├── Model:       ", sim.model)
    println(io, "├── Flux:        ", sim.flux)
    println(io, "├── Mode:        ", sim.mode)
    println(io, "├── Savedir:     ", sim.savedir)
    show_time_fields(io, sim.time)
    println(io, "├── Cache:       ", sim.cache_initialized ? "initialized" : "not initialized")
    print(io,   "└── Save flux:   ", sim.save_input_flux)
end

function show_time_fields(io::IO, time::RefinedTimeGrid)
    println(io, "├── duration:    ", time.duration, " s")
    println(io, "├── dt request:  ", time.dt_requested, " s")
    println(io, "├── dt resolved: ", time.dt_resolved, " s")
    println(io, "├── CFL factor:  ", time.CFL_factor)
    println(io, "├── n_loop:      ", time.n_loop)
end
function show_time_fields(io::IO, time::UniformTimeGrid)
    println(io, "├── duration:    ", time.duration, " s")
    println(io, "├── dt:          ", time.dt, " s")
    println(io, "├── n_steps:     ", time.n_steps)
end
function show_time_fields(io::IO, ::SingleStepConfig)
    println(io, "├── Time:        single-step")
end
