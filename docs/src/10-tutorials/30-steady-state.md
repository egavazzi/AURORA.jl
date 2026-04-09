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
    180:-10:0,     # pitch-angle bin edges [°]
    3000,          # maximum energy [eV]
    msis_file,
    iri_file,
    13             # magnetic field angle to zenith [°]
)
```

## Flat spectrum with energy cutoff

A [`FlatSpectrum`](@ref) produces a uniform differential number flux above a given minimum
energy. This is useful for studying the transport of quasi-monoenergetic beams:

```@example steady_state
flux = InputFlux(
    FlatSpectrum(1e-3; E_min=2900);    # 1 mW/m², only above 2900 eV
    beams=1:2                          # two most field-aligned downward beams
)

savedir = mkpath(joinpath("data", "steady_state_example"))
sim = AuroraSimulation(model, flux, savedir)
run!(sim)
nothing # hide
```

## Analysis

After the simulation completes, compute derived quantities:

```@example steady_state
# Volume excitation rates for optical emissions (e.g. N₂⁺ 1NG, O 557.7 nm, ...)
make_volume_excitation_file(sim)

# Precipitated current density
make_current_file(sim)

# Input flux boundary condition (for verification)
make_Ie_top_file(sim)
```

```@example steady_state
readdir(savedir)
```

## Gaussian and Maxwellian spectra

For more realistic energy distributions, use [`GaussianSpectrum`](@ref) or
[`MaxwellianSpectrum`](@ref):

```@example steady_state
# Gaussian: 10 mW/m², peaked at 2 keV with 300 eV width
flux_gauss = InputFlux(
    GaussianSpectrum(1e-2, 2000.0, 300.0);
    beams=1:2
)
```

```@example steady_state
# Maxwellian: 10 mW/m², characteristic energy 1 keV
flux_maxw = InputFlux(
    MaxwellianSpectrum(1e-2, 1000.0);
    beams=1:2
)
```

See the [Input flux](@ref "Input Flux") tutorial for the full list of
spectrum and modulation types.

For the general time-dependent case, including flickering input fluxes, see
[Time-Dependent Simulation](@ref "Time-Dependent Simulation").
