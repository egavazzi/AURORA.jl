using Dates: now, Dates

"""
    AuroraOutputManager(savedir; overwrite=false, compress=true)

Configuration struct that controls where and how simulation output is written.

# Arguments
- `savedir`: directory where output files will be written (created by `run!` if absent).
  May be a relative or absolute path. If empty (or only whitespace), output goes to
  `backup/<yyyymmdd-HHMM>/` in the current working directory.

# Keyword arguments
- `overwrite::Bool=false`: if `false` (the default), `run!()` errors when
  `simulation_data.nc` already exists in `savedir`; set to `true` to allow overwriting.
- `compress::Bool=true`: enable zlib compression (`deflatelevel=4`) on all NetCDF variables.

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
out = AuroraOutputManager("my_run"; compress=false)
sim = AuroraSimulation(model, flux, out; mode=TimeDependentMode(duration=0.5, dt=0.001))

# Convenience — pass a plain String and defaults apply:
sim = AuroraSimulation(model, flux, "my_run"; mode=SteadyStateMode())
```
"""
struct AuroraOutputManager
    savedir::String
    overwrite::Bool
    compress::Bool
end

AuroraOutputManager(savedir; overwrite=false, compress=true) =
    AuroraOutputManager(resolve_savedir(savedir), overwrite, compress)

"""
    resolve_savedir(savedir) -> String

Normalise a user-supplied `savedir`. A non-empty path (relative or absolute) is returned
unchanged; an empty or whitespace-only path falls back to `backup/<yyyymmdd-HHMM>` in the
current working directory.
"""
function resolve_savedir(savedir)
    s = String(savedir)
    if isempty(s) || !occursin(r"[^ ]", s)
        return joinpath("backup", Dates.format(now(), "yyyymmdd-HHMM"))
    end
    return s
end

function Base.show(io::IO, out::AuroraOutputManager)
    print(io, "AuroraOutputManager(\"", out.savedir, "\")")
end

function Base.show(io::IO, ::MIME"text/plain", out::AuroraOutputManager)
    println(io, "AuroraOutputManager:")
    println(io, "├── savedir:    ", out.savedir)
    println(io, "├── overwrite:  ", out.overwrite)
    print(io,   "└── compress:   ", out.compress)
end
