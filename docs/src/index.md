```@meta
CurrentModule = AURORA
```

# AURORA.jl

AURORA.jl is a time-dependent multi-stream electron transport code.

It is suitable for modeling the transport of ionospheric electron fluxes when 
the electron precipitation varies rapidly, typically on sub-second timescales.

## Quick start

```@example quick-start
using AURORA

# Set up the physical model
msis_file = find_msis_file()
iri_file  = find_iri_file()
model = AuroraModel([100, 600], 180:-10:0, 1000, msis_file, iri_file, 13)

# Define a flickering input flux
flux = InputFlux(
    FlatSpectrum(1e-2; E_min=100),
    SinusoidalFlickering(5.0);
    beams=1,
    z_source=3000.0,
)

# Create and run a short time-dependent simulation
savedir = mkpath(joinpath("data", "my_first_simulation"))
sim = AuroraSimulation(model, flux, 0.2, 0.01, savedir; CFL_number=128)
run!(sim)
```

See the [Time-Dependent Simulation](@ref) tutorial for a full walkthrough.

## Documentation outline

- **[Tutorials](@ref "Installation")** — step-by-step guides: installation, first simulation,
  steady-state and time-dependent examples, custom input flux.
- **[Physics](@ref "Transport Equation")** — the equations and physical models behind AURORA.
- **[Implementation](@ref "Architecture")** — internal architecture, call chain, matrix
    assembly, and solver details.
- **[API Reference](@ref "Model")** — docstrings for all exported types and functions.

## Citing AURORA

If you use AURORA.jl in your work, please cite it as follows. Note that
`10.5281/zenodo.11238620` refers to the general project and all its versions.
Version-specific DOIs can be found on Zenodo.

```bibtex
@software{Gavazzi_AURORA_jl_2026,
    author = {Gavazzi, Etienne and Gustavsson, Björn},
    doi = {10.5281/zenodo.11238620},
    license = {GPL-3.0},
    month = mar,
    title = {{AURORA.jl}},
    url = {https://github.com/egavazzi/AURORA.jl},
    version = {0.7.0},
    year = {2026}
}
```

