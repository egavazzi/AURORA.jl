# Cross Sections

!!! info "Under construction"
    This section is a placeholder. Detailed tables and plots of the cross-section data
    will be added in a future update.

AURORA includes electron-impact cross sections for three neutral species in the
upper atmosphere:

| Species | Elastic | Inelastic (excitation) | Ionization |
|---------|---------|----------------------|------------|
| N₂ | ✓ | Multiple rotational, vibrational, and electronic states | ✓ |
| O₂ | ✓ | Multiple states | ✓ |
| O  | ✓ | Multiple states | ✓ |

## Data sources

The cross-section data are primarily based on:

- **Itikawa, Y.** (2006). Cross sections for electron collisions with nitrogen molecules.
  *J. Phys. Chem. Ref. Data*, 35(1), 31–53.
- **Itikawa, Y.** (2009). Cross sections for electron collisions with oxygen molecules.
  *J. Phys. Chem. Ref. Data*, 38(1), 1–20.
- **Itikawa, Y. & Ichimura, A.** (1990). Cross sections for collisions of electrons and
  photons with atomic oxygen. *J. Phys. Chem. Ref. Data*, 19(3), 637–651.

## Current scope in AURORA

The current model includes cross sections for N₂, O₂, and O, together with selected
optical-emission cross sections used to compute excitation rates. These data are loaded
internally when an [`AuroraModel`](@ref) is constructed.

## Emission cross sections

AURORA also includes cross sections for specific auroral optical emissions, used by
[`make_volume_excitation_file`](@ref) to compute volume excitation rates. These include
transitions such as:

- N₂⁺ first negative band (1NG, 3914 Å)
- (4278 Å)
- O(¹S) green line (5577 Å)
- O(¹D) red line (6300 Å)
- 7774 Å
- 8446 Å
- N₂ first positive band (1PG)
- N₂ second positive band (2PG)
