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

    n_intervals = Float64(duration) / Float64(dt)
    n_intervals_round = round(Int, n_intervals)
    isapprox(n_intervals, n_intervals_round; atol=1e-10, rtol=1e-10) ||
        error("duration ($duration) must be an integer multiple of dt ($dt)")

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

struct RefinedTimeGrid{T<:AbstractRange{Float64}} <: AbstractTimeConfig
    duration::Float64
    dt_requested::Float64
    dt_resolved::Float64
    CFL_factor::Int
    t::T
    n_loop::Int
    n_t_per_loop::Int
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

    # Calculate timesteps per loop and actual resolved timestep (can be different from the saving dt)
    n_t_per_loop = (length(t_resolved) - 1) ÷ n_loop_resolved + 1
    dt_resolved = length(t_resolved) > 1 ? Float64(t_resolved[2] - t_resolved[1]) : Float64(dt)

    # Return fully configured time grid
    return RefinedTimeGrid(Float64(duration), Float64(dt), dt_resolved,
                            CFL_factor, t_resolved, n_loop_resolved, n_t_per_loop)
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
    print(io,   "└── steps / loop:  ", time.n_t_per_loop)
end


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
mutable struct AuroraSimulation{M<:AuroraModel, F<:InputFlux, S<:AbstractMode}
    const model::M
    const flux::F
    const mode::S
    const savedir::String
    const time::AbstractTimeConfig
    const save_input_flux::Bool
    cache::Union{Nothing, SimulationCache}
end

function AuroraSimulation(model::AuroraModel, flux::InputFlux, savedir;
                          mode::AbstractMode=SteadyStateMode(),
                          save_input_flux=true)
    time = _build_time_config(model, mode)
    return AuroraSimulation(model, flux, mode, String(savedir), time, save_input_flux, nothing)
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
    _show_time_fields(io, sim.time)
    println(io, "├── Cache:       ", sim.cache === nothing ? "not initialized" : "initialized")
    print(io,   "└── Save flux:   ", sim.save_input_flux)
end

function _show_time_fields(io::IO, time::RefinedTimeGrid)
    println(io, "├── duration:    ", time.duration, " s")
    println(io, "├── dt request:  ", time.dt_requested, " s")
    println(io, "├── dt resolved: ", time.dt_resolved, " s")
    println(io, "├── CFL factor:  ", time.CFL_factor)
    println(io, "├── n_loop:      ", time.n_loop)
end
function _show_time_fields(io::IO, time::UniformTimeGrid)
    println(io, "├── duration:    ", time.duration, " s")
    println(io, "├── dt:          ", time.dt, " s")
    println(io, "├── n_steps:     ", time.n_steps)
end
function _show_time_fields(io::IO, ::SingleStepConfig)
    println(io, "├── Time:        single-step")
end
