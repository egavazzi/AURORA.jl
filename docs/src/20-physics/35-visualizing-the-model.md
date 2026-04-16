# Visualizing the Model

After building an [`AuroraModel`](@ref), [`plot_model`](@ref) can be called to produce a set of diagnostic
figures that let you inspect the atmosphere, energy grid, cross-sections, phase functions,
and scattering data before running a simulation.

## Quick start

```@example plot_model
using AURORA
using CairoMakie
CairoMakie.activate!() # hide

altitude_lims = [100, 600]
θ_lims = 180:-10:0
E_max = 3000

msis_file = find_msis_file(year=2005, month=10, day=8, hour=22, minute=0,
                           lat=70, lon=19, height=85:1:700)
iri_file = find_iri_file(year=2005, month=10, day=8, hour=22, minute=0,
                         lat=70, lon=19, height=85:1:700)

model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file)
```

Generate all figures at once:

```@example plot_model
figs = plot_model(model)
```

Or select specific panels:

```@example plot_model
figs_subset = plot_model(model; panels=[:atmosphere, :cross_sections])
```

## Panels

### Atmosphere

Neutral density profiles (N₂, O₂, O) from MSIS and the electron density profile from IRI,
both interpolated onto the model altitude grid.

```@example plot_model
figs[:atmosphere]
```

### Energy levels

Excitation thresholds (vibrational, rotational, electronic, and ionization) for each
neutral species.

```@example plot_model
figs[:energy_levels]
```

### Energy grid

Variation of the energy bin width with energy. The grid is finer at low energies and
coarser at high energies.

```@example plot_model
figs[:energy_grid]
```

### Cross-sections

Electron-impact cross-sections for all excited/ionized states of N₂, O₂, and O.
Solid lines are excitation states and dashed lines are ionization states.
Rotational/vibrational excitations are plotted with thinner lines.

```@example plot_model
figs[:cross_sections]
```

### Phase functions

The scattering phase functions (elastic and inelastic) show, for each species, the probability distribution 
for scattering angles as a function of energy.

```@example plot_model
figs[:phase_functions]
```

### Scattering matrices

Pitch-angle transition probability matrices for each pitch-angle stream. Each panel shows
the probability that an electron with an initial pitch-angle "from-angle" will end up in the 
stream/beam corresponding to the panel, after scattering with an angle "scattering-angle".

```@example plot_model
figs[:scattering]
```
