# Energy Degradation

!!! info "Under construction"
    This section is a placeholder. A detailed description of the cascading and
    scattering mechanisms will be added in a future update.

When a precipitating electron collides with a neutral atom or molecule, it can lose
energy through several processes. AURORA tracks these energy transfers to build the
source term ``Q`` that feeds lower-energy bins.

## Processes

### Inelastic scattering

An electron excites a neutral to a higher electronic, vibrational, or rotational state,
losing a discrete amount of energy equal to the excitation threshold. The scattered
electron continues at reduced energy with a new pitch angle determined by the
differential cross section.

### Ionization

An electron ionizes a neutral, producing:

1. A **secondary electron** — ejected from the neutral.
2. A **degraded primary** — the incident electron, which continues with energy reduced by
   the ionization potential plus the energy given to the secondary electron.

Both the secondary and degraded primary electrons contribute to the source term for
lower energy bins.

### Electron-electron collisions

Energetic electrons also lose energy through Coulomb collisions with the thermal electron
population. This continuous energy loss is modelled using the Swartz & Nisbet (1971)
formula and transfers energy from bin ``iE`` to bin ``iE-1``.

## Cascading transfer matrices

In AURORA, the redistribution from ionization is represented using transfer matrices for
each neutral species (N₂, O₂, O). These matrices map the ionization rate at a high energy
``E'`` to the production of secondary and degraded-primary electrons at lower energies.

See [Internals](@ref) for how these couplings are assembled in the numerical solver.
