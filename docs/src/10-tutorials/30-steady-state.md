# [Steady-State Simulation](@id Steady-State-Simulation)

AURORA can also run steady-state simulations. In this mode, it computes the final 
equilibrium electron distribution for a constant (time-independent) precipitating
input. This is useful when you only care about the long-term response and don't 
need to resolve the precise temporal dynamics.

## Setup

```@example steady_state
using AURORA

# Atmospheric model — default VISIONS-2 conditions
msis_file = find_msis_file()
iri_file  = find_iri_file()

model = AuroraModel(
    [100, 600],    # altitude limits [km]
    180:-15:0,     # pitch-angle bin edges [°]
    3000,          # maximum energy [eV]
    msis_file,
    iri_file,
    13             # magnetic field angle to zenith [°]
)
```

## Maxwellian input spectrum

A [`MaxwellianSpectrum`](@ref) provides a smooth, physically motivated energy distribution
set by a characteristic energy (formula from the Meier et al. 1989 paper):

```@example steady_state
flux = InputFlux(
    MaxwellianSpectrum(1e-3, 1000.0);  # 1 mW/m², characteristic energy 1 keV
    beams=1:2                          # two most field-aligned downward beams
)
```

## Create and run the simulation

```@example steady_state
savedir = mkpath(joinpath("data", "steady_state_example"))

sim = AuroraSimulation(model, flux, savedir)
```

```@setup steady_state
# Initialize the simulation to suppress the verbose "Load/Calculate..." messages
# from appearing in the example output.
initialize!(sim)
```

```@example steady_state
run!(sim)
nothing # hide
```

## Post-process

As for [time-dependent](@ref "Time-Dependent Simulation") simulations, it is possible to compute derived quantities from the raw
electron flux and save them alongside the simulation output:

```@example steady_state
make_Ie_top_file(sim)              # boundary condition (input flux applied at top)
make_volume_excitation_file(sim)   # volumetric excitation rates for optical emissions
make_column_excitation_file(sim)   # column-integrated excitation rates (steady-state scalar)
make_current_file(sim)             # field-aligned electron currents and energy fluxes
make_heating_rate_file(sim)        # electron heating rates
nothing # hide
```

```@example steady_state
readdir(savedir)
```

## Visualize

AURORA.jl provides helper plotting functions through a [Makie](https://github.com/MakieOrg/Makie.jl) extension.
Install and load a Makie backend to access them (more information in [`Visualization`](@ref "Visualization")).

### Input flux

[`plot_input`](@ref) can be used to inspect the prescribed input spectrum directly from the simulation object:

```@example steady_state
using CairoMakie

fig = plot_input(sim)
```

### Volume excitation profile

[`plot_excitation!`](@ref) can be used to visualize the altitude profile of the volume excitation rate for a given emission line. 

```@example steady_state
using CairoMakie
vol = load_volume_excitation(savedir)

fig = Figure()
ax  = Axis(fig[1, 1];
           xlabel = "Volume emission rate (photons/m³/s)",
           ylabel = "Altitude (km)",
           xscale = log10,
           title = "4278 Å")
plot_excitation!(ax, vol; field = :Q4278)
fig
```

## Next steps

- [Time-Dependent Simulation](@ref "Time-Dependent Simulation") — run the general time-dependent case.
- [Input flux](@ref "Input Flux") — explore all spectrum and modulation types.
- [Post-processing & analysis](@ref "Post-Processing") — detailed walkthrough of all analysis functions.
