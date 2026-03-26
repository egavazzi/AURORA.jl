struct ResolvedTimeGrid{T<:AbstractRange{Float64}}
    t_total::Float64
    dt_requested::Float64
    dt_resolved::Float64
    CFL_number::Float64
    CFL_factor::Int
    t::T
    n_loop::Int
    max_memory_gb::Float64
    n_t_per_loop::Int
end

function ResolvedTimeGrid(model::AuroraModel, t_total, dt;
                          CFL_number=64, n_loop=nothing, max_memory_gb=8.0)
    t_total > 0 || error("t_total must be positive, got $t_total")
    dt > 0 || error("dt must be positive, got $dt")
    dt <= t_total || error("dt ($dt) must be ≤ t_total ($t_total)")
    CFL_number > 0 || error("CFL_number must be positive, got $CFL_number")
    max_memory_gb > 0 || error("max_memory_gb must be positive, got $max_memory_gb")
    if !isnothing(n_loop)
        n_loop > 0 || error("n_loop must be positive, got $n_loop")
    end

    z = model.altitude_grid.h
    n_μ = length(model.pitch_angle_grid.μ_center)
    n_E = model.energy_grid.n
    v_max = v_of_E(maximum(model.energy_grid.E_centers))

    t_resolved, CFL_factor = CFL_criteria(t_total, dt, z, v_max, CFL_number)
    n_loop_resolved = isnothing(n_loop) ? calculate_n_loop(t_resolved, length(z), n_μ, n_E;
                                                           max_memory_gb=max_memory_gb) : Int(n_loop)
    check_n_loop(n_loop_resolved, length(z), n_μ, length(t_resolved), n_E)

    n_t_per_loop = (length(t_resolved) - 1) ÷ n_loop_resolved + 1
    dt_resolved = length(t_resolved) > 1 ? Float64(t_resolved[2] - t_resolved[1]) : Float64(dt)

    return ResolvedTimeGrid(Float64(t_total), Float64(dt), dt_resolved, Float64(CFL_number),
                            CFL_factor, t_resolved, n_loop_resolved,
                            Float64(max_memory_gb), n_t_per_loop)
end

function Base.show(io::IO, time::ResolvedTimeGrid)
    print(io, "ResolvedTimeGrid(t_total=$(time.t_total), dt=$(time.dt_requested), n_loop=$(time.n_loop))")
end

function Base.show(io::IO, ::MIME"text/plain", time::ResolvedTimeGrid)
    println(io, "ResolvedTimeGrid:")
    println(io, "├── Total time:    ", time.t_total, " s")
    println(io, "├── Requested dt:  ", time.dt_requested, " s")
    println(io, "├── Resolved dt:   ", time.dt_resolved, " s")
    println(io, "├── CFL number:    ", time.CFL_number)
    println(io, "├── CFL factor:    ", time.CFL_factor)
    println(io, "├── n_loop:        ", time.n_loop)
    println(io, "├── steps per loop:    ", time.n_t_per_loop)
    print(io,   "└── memory budget: ", time.max_memory_gb, " GB")
end

"""
    AuroraSimulation{M<:AuroraModel, F<:InputFlux}

A complete simulation configuration, bundling the physical model, input flux,
output directory, resolved time grid, and lazy simulation cache.

Use [`initialize!`](@ref) to allocate the workspace explicitly, or call [`run!`](@ref)
directly to auto-initialize and execute the simulation.
"""
mutable struct AuroraSimulation{M<:AuroraModel, F<:InputFlux}
    const model::M
    const flux::F
    const savedir::String
    const time::Union{ResolvedTimeGrid, Nothing}
    const save_input_flux::Bool
    cache::Union{Nothing, TransportCache}
end

function AuroraSimulation(model::AuroraModel, flux::InputFlux, t_total, dt, savedir;
                          CFL_number=64, n_loop=nothing, max_memory_gb=8.0,
                          save_input_flux=true)
    time = ResolvedTimeGrid(model, t_total, dt; CFL_number, n_loop, max_memory_gb)
    return AuroraSimulation(model, flux, String(savedir), time, save_input_flux, nothing)
end

function AuroraSimulation(model::AuroraModel, flux::InputFlux, savedir;
                          save_input_flux=true)
    return AuroraSimulation(model, flux, String(savedir), nothing, save_input_flux, nothing)
end

function Base.show(io::IO, sim::AuroraSimulation)
    mode = isnothing(sim.time) ? "steady-state" : "time-dependent"
    print(io, "AuroraSimulation($mode)")
end

function Base.show(io::IO, ::MIME"text/plain", sim::AuroraSimulation)
    mode = isnothing(sim.time) ? "Steady-state" : "Time-dependent"
    println(io, "AuroraSimulation ($mode):")
    println(io, "├── Model:       ", sim.model)
    println(io, "├── Flux:        ", sim.flux)
    println(io, "├── Savedir:     ", sim.savedir)
    if isnothing(sim.time)
        println(io, "├── Time:        steady-state")
    else
        println(io, "├── t_total:     ", sim.time.t_total, " s")
        println(io, "├── dt request:  ", sim.time.dt_requested, " s")
        println(io, "├── dt resolved: ", sim.time.dt_resolved, " s")
        println(io, "├── CFL:         ", sim.time.CFL_number)
        println(io, "├── n_loop:      ", sim.time.n_loop)
    end
    println(io, "├── Cache:       ", sim.cache === nothing ? "not initialized" : "initialized")
    print(io,   "└── Save flux:   ", sim.save_input_flux)
end


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
    step!(sim)
    finalize!(sim)
    return sim
end
