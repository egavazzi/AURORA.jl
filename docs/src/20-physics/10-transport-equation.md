# Transport Equation

!!! info "Under construction"
    This section is a placeholder. A detailed derivation and discussion of the multi-stream
    transport equation will be added in a future update.

AURORA solves the time-dependent, multi-stream electron transport equation in the
Earth's ionosphere. The equation describes how precipitating auroral electrons propagate
along the magnetic field, scatter off neutral atoms and molecules, lose energy through
inelastic collisions, and produce secondary electrons through ionization.

## The equation

```math
\frac{\partial I_e}{\partial t} =
    -\mu \, v \frac{\partial I_e}{\partial z}
    + \frac{\partial}{\partial z}\left(D \frac{\partial I_e}{\partial z}\right)
    - A \cdot I_e
    + B \cdot I_e
    + Q
```

| Term | Physical meaning |
|------|-----------------|
| ``-\mu v \, \partial I_e / \partial z`` | Field-aligned streaming |
| ``\partial_z(D \, \partial_z I_e)`` | Spatial diffusion (pitch-angle mixing) |
| ``-A \cdot I_e`` | Losses (collisions with neutrals + thermal electrons) |
| ``+B \cdot I_e`` | Pitch-angle scattering (elastic + inelastic) |
| ``+Q`` | Sources from energy cascading (ionization secondaries) |

The independent variables are altitude ``z``, pitch-angle cosine ``\mu``, time ``t``, and
energy ``E``. The loss and source terms couple the solution across the energy grid because
higher-energy electrons can degrade and produce lower-energy secondaries.

For how this equation is discretized and solved in AURORA, see [Internals](@ref).

## References

- Gustavsson, B. (2022). *Time-Dependent Electron Transport I: Modelling of Supra-Thermal
  Electron Bursts Modulated at 5–10 Hz With Implications for Flickering Aurora*. Journal of
  Geophysical Research: Space Physics, 127(6), e2019JA027608.
- Gavazzi, E. (2022). *The effects of time-variation of electron fluxes from the auroral
  ionosphere on M-I coupling* [Master thesis, UiT Norges arktiske universitet].
