"""
    AuroraOutputManager(savedir; overwrite=false, compress=true, save_input_flux=true)

Configuration struct that controls where and how simulation output is written.

# Arguments
- `savedir::String`: directory where output files will be written (created if absent).

# Keyword arguments
- `overwrite::Bool=false`: if `false` (the default), `run!()` errors when
  `simulation_data.nc` already exists in `savedir`; set to `true` to allow overwriting.
- `compress::Bool=true`: enable zlib compression (`deflatelevel=4`) on all NetCDF variables.
- `save_input_flux::Bool=true`: write the top-of-atmosphere input flux (precipitation) as a
  separate `Ie_input` variable (on its own `time_input` axis) in `simulation_data.nc`.

# Output layout
```
savedir/
├── config.toml               # scalar simulation parameters
├── inputs/
│   ├── atmosphere.nc         # neutral/ionosphere density profiles
│   └── physics_state.jld2   # full AuroraModel (species + physics matrices)
├── simulation_data.nc        # time-dependent flux output (appended per loop)
└── analysis/                 # derived quantities (written by make_* functions)
```

# Examples
```julia
# Full control:
out = AuroraOutputManager("my_run"; compress=false, save_input_flux=false)
sim = AuroraSimulation(model, flux, out; mode=TimeDependentMode(duration=0.5, dt=0.001))

# Convenience — pass a plain String and defaults apply:
sim = AuroraSimulation(model, flux, "my_run"; mode=SteadyStateMode())
```
"""
struct AuroraOutputManager
    savedir::String
    overwrite::Bool
    compress::Bool
    save_input_flux::Bool
end

AuroraOutputManager(savedir; overwrite=false, compress=true, save_input_flux=true) =
    AuroraOutputManager(String(savedir), overwrite, compress, save_input_flux)

function Base.show(io::IO, out::AuroraOutputManager)
    print(io, "AuroraOutputManager(\"", out.savedir, "\")")
end

function Base.show(io::IO, ::MIME"text/plain", out::AuroraOutputManager)
    println(io, "AuroraOutputManager:")
    println(io, "├── savedir:          ", out.savedir)
    println(io, "├── overwrite:        ", out.overwrite)
    println(io, "├── compress:         ", out.compress)
    print(io,   "└── save_input_flux:  ", out.save_input_flux)
end
