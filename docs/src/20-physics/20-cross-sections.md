# Cross Sections

!!! info "WIP"
    This section is under construction.

AURORA includes electron-impact cross sections for the three main neutral species
in the upper atmosphere: **N₂**, **O₂**, and **atomic O**. For each species, elastic and
inelastic (excitation of rotational, vibrational, electronic states and ionization) 
cross-sections are included. 

Each species carries its inelastic collisions as a vector of [`CollisionChannel`](@ref)s —
a name, a cross-section function σ(E) in m², the energy loss in eV, the number of secondary
electrons ejected (0, 1 or 2), and a free-form `source` string recording where the data come
from. Its elastic cross section and its secondary-electron energy distribution sit alongside,
in the `elastic_cross_section` and `secondary_law` fields of [`NeutralSpecies`](@ref).

At `initialize!(model)` the channel table is evaluated on the model energy grid into the
`cross_sections` matrix `[n_levels × n_E]` and the `excitation_levels` matrix
`[n_levels × 2]`, and the ionizing channels define the cascading thresholds. Row 1 of both
matrices is the elastic channel and row `i + 1` is `channels[i]`, so the rows cannot fall out
of step and the cascading thresholds cannot disagree with the energy losses charged to the
primary electron.

The built-in tables live in `src/physics/cross_sections/channels_N2.jl`, `channels_O2.jl` and
`channels_O.jl`, and are returned by [`default_channels`](@ref),
[`default_elastic_cross_section`](@ref) and [`default_secondary_law`](@ref). A run writes its
tables to `inputs/collision_channels.toml` for inspection.

## Data sources

The cross-section data are primarily based on:

- **Itikawa, Y.** (2006). Cross sections for electron collisions with nitrogen molecules.
  *J. Phys. Chem. Ref. Data*, 35(1), 31–53.
- **Itikawa, Y.** (2009). Cross sections for electron collisions with oxygen molecules.
  *J. Phys. Chem. Ref. Data*, 38(1), 1–20.
- **Itikawa, Y. & Ichimura, A.** (1990). Cross sections for collisions of electrons and
  photons with atomic oxygen. *J. Phys. Chem. Ref. Data*, 19(3), 637–651.

## Emission cross sections

In addition to the collision cross sections used by the transport solver, AURORA includes
cross sections for specific auroral optical emissions. These are used by
[`make_volume_excitation_file`](@ref) to compute volume excitation rates. The following
emission lines are currently implemented:

| Emission | Source | Wavelength | Function |
|----------|--------|------------|----------|
| N₂⁺ 1NG (0-1) | e + N₂ | 4278 Å | `excitation_4278` |
| N₂ 1PG (4–1, 5–2) | e + N₂ | 6730 Å | `excitation_6730_N2` |
| OI 7774 Å | e + O | 7774 Å | `excitation_7774_O` |
| OI 7774 Å | e + O₂ | 7774 Å | `excitation_7774_O2` |
| OI 8446 Å | e + O | 8446 Å | `excitation_8446_O` |
| OI 8446 Å | e + O₂ | 8446 Å | `excitation_8446_O2` |
| O(¹S) green line | e + O | 5577 Å | `excitation_O1S` |
| O(¹D) red line | e + O | 6300 Å | `excitation_O1D` |
