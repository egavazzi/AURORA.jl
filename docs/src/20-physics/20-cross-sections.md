# Cross Sections

!!! info "WIP"
    This section is under construction.

AURORA includes electron-impact cross sections for the three main neutral species
in the upper atmosphere: **N₂**, **O₂**, and **atomic O**. For each species, elastic and
inelastic (excitation of rotational, vibrational, electronic states and ionization) 
cross-sections are included. 

These data are loaded internally when an [`AuroraModel`](@ref) is constructed.

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
| O(¹S) green line | e + O | 5577 Å | `excitation_O1S` |
| O(¹D) red line | e + O | 6300 Å | `excitation_O1D` |
| N₂ 1PG (4–1, 5–2) | e + N₂ | 6730 Å | `excitation_6730_N2` |
| OI 7774 Å | e + O | 7774 Å | `excitation_7774_O` |
| OI 7774 Å | e + O₂ | 7774 Å | `excitation_7774_O2` |
| OI 8446 Å | e + O | 8446 Å | `excitation_8446_O` |
| OI 8446 Å | e + O₂ | 8446 Å | `excitation_8446_O2` |
