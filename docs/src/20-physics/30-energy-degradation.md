# Energy Degradation

!!! info "WIP"
    This section is under construction.

When a precipitating electron collides with a neutral atom or molecule, it can lose
energy through several processes. AURORA tracks these energy transfers to build the
source term ``Q`` that feeds lower-energy bins.

## Processes

### Electron-electron collisions

Energetic electrons also lose energy through Coulomb collisions with the thermal electron
population. This continuous energy loss is modelled using the Swartz & Nisbet (1971)
formula and transfers energy from bin ``iE`` to bin ``iE-1``.

```math
L_e(E, z) = \frac{3.0271 \times 10^{-10} \, n_e(z)^{0.97}}{E^{0.44}} \left(\frac{E - E_e(z)}{E - 0.53\,E_e(z)}\right)^{2.36} \frac{1}{v(E)}
```

with thermal electron energy

```math
E_e(z) = \frac{k_B}{q_e} T_e(z)
```

where ``L_e`` is the suprathermal-electron energy loss rate [eV m⁻¹], ``n_e`` is the
ambient electron density [m⁻³], ``T_e`` is the thermal electron temperature [K], and
``v(E)`` is the speed of an electron with energy ``E``.


### Inelastic scattering

An electron excites a neutral to a higher electronic, vibrational, or rotational state,
losing a discrete amount of energy equal to the excitation threshold. The scattered
electron continues at reduced energy with a new pitch angle determined by the
differential cross section.

### Ionization

An electron ionizes a neutral, producing:

1. A **secondary electron** — ejected from the neutral. Secondary electrons are emitted
   **isotropically** over all pitch angles, i.e. they are redistributed uniformly across
   all beams weighted by their solid angle. Their energy distribution is
   species-dependent and evaluated from a species-dependent analytic law before being
   rebinned onto the simulation energy grid.
2. A **degraded primary** — the incident electron, which continues at reduced energy
   (reduced by the ionization potential plus the energy given to the secondary) and in the
   **same pitch-angle beam** as the original electron. Its lower-energy redistribution is
   obtained from the pre-computed cascading transfer matrices.

Both the secondary and degraded primary electrons contribute to the source term for
lower energy bins.
