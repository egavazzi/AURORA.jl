"""
    AuroraSimulation{M<:AuroraModel, F<:InputFlux}

A complete simulation configuration, bundling the physical model, input flux,
output directory, and solver settings.

Use [`run!`](@ref) to execute the simulation.

# Fields
- `model`: the [`AuroraModel`](@ref) describing the grids and atmosphere/ionosphere
- `flux`: the [`InputFlux`](@ref) describing the incoming electron precipitation
- `savedir`: path to the directory where results will be saved
- `t_total`: total simulation time (s), or `nothing` for steady-state
- `dt`: output time step (s), or `nothing` for steady-state
- `CFL_number`: CFL multiplier for time sub-sampling (default: `64`)
- `n_loop`: number of time-chunks, or `nothing` for automatic
- `max_memory_gb`: memory budget (GB) for automatic `n_loop` computation
- `save_input_flux`: whether to save the computed top-boundary flux before the simulation

# Constructors
```julia
# Time-dependent simulation
AuroraSimulation(model, flux, t_total, dt, savedir;
                 CFL_number=64, n_loop=nothing, max_memory_gb=8, save_input_flux=true)

# Steady-state simulation
AuroraSimulation(model, flux, savedir; save_input_flux=true)
```

# Examples
```julia
model = AuroraModel([100, 600], 180:-10:0, 1000, msis_file, iri_file, 13)
flux  = InputFlux(FlatSpectrum(1e-2; E_min=100), SinusoidalFlickering(5.0);
                  beams=1, z_source=3000.0)
savedir = make_savedir("backup", "my_experiment")

# Time-dependent
sim = AuroraSimulation(model, flux, 0.5, 0.001, savedir; CFL_number=128)
run!(sim)

# Steady-state
sim = AuroraSimulation(model, flux_ss, savedir)
run!(sim)
```
"""
struct AuroraSimulation{M<:AuroraModel, F<:InputFlux}
    model::M
    flux::F
    savedir::String
    t_total::Union{Float64, Nothing}
    dt::Union{Float64, Nothing}
    CFL_number::Float64
    n_loop::Union{Int, Nothing}
    max_memory_gb::Float64
    save_input_flux::Bool
end

# Time-dependent constructor
function AuroraSimulation(model::AuroraModel, flux::InputFlux, t_total, dt, savedir;
                          CFL_number=64, n_loop=nothing, max_memory_gb=8.0,
                          save_input_flux=true)
    t_total > 0 || error("t_total must be positive, got $t_total")
    dt > 0 || error("dt must be positive, got $dt")
    dt <= t_total || error("dt ($dt) must be ≤ t_total ($t_total)")
    CFL_number > 0 || error("CFL_number must be positive, got $CFL_number")
    max_memory_gb > 0 || error("max_memory_gb must be positive, got $max_memory_gb")
    if !isnothing(n_loop)
        n_loop > 0 || error("n_loop must be positive, got $n_loop")
    end
    AuroraSimulation(model, flux, String(savedir), Float64(t_total), Float64(dt),
                     Float64(CFL_number), n_loop, Float64(max_memory_gb), save_input_flux)
end

# Steady-state constructor
function AuroraSimulation(model::AuroraModel, flux::InputFlux, savedir;
                          save_input_flux=true)
    AuroraSimulation(model, flux, String(savedir), nothing, nothing,
                     0.0, nothing, 0.0, save_input_flux)
end

function Base.show(io::IO, sim::AuroraSimulation)
    mode = isnothing(sim.t_total) ? "steady-state" : "time-dependent"
    print(io, "AuroraSimulation($mode)")
end

function Base.show(io::IO, ::MIME"text/plain", sim::AuroraSimulation)
    mode = isnothing(sim.t_total) ? "Steady-state" : "Time-dependent"
    println(io, "AuroraSimulation ($mode):")
    println(io, "├── Model:    ", sim.model)
    println(io, "├── Flux:     ", sim.flux)
    println(io, "├── Savedir:  ", sim.savedir)
    if !isnothing(sim.t_total)
        println(io, "├── t_total:  ", sim.t_total, " s")
        println(io, "├── dt:       ", sim.dt, " s")
        println(io, "├── CFL:      ", sim.CFL_number)
        n_str = isnothing(sim.n_loop) ? "auto (≤ $(sim.max_memory_gb) GB)" : string(sim.n_loop)
        println(io, "├── n_loop:   ", n_str)
    end
    print(io,   "└── Save flux: ", sim.save_input_flux)
end


"""
    run!(sim::AuroraSimulation)

Execute the simulation described by `sim`.

Dispatches to time-dependent or steady-state solver based on whether `t_total` is set.

# Examples
```julia
sim = AuroraSimulation(model, flux, 0.5, 0.001, savedir)
run!(sim)
```
"""
function run!(sim::AuroraSimulation)
    if isnothing(sim.t_total)
        calculate_e_transport_steady_state(sim)
    else
        calculate_e_transport(sim)
    end
end
