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
