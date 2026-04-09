# [Analysis](@id Analysis)

Post-processing functions that compute derived quantities from raw simulation output.
All functions write their results to the simulation's save directory.

!!! note "Compatibility"

    | Function | Steady-State results | Time-Dependent results |
    |----------|:------------:|:--------------:|
    | `make_volume_excitation_file` | ✓ | ✓ |
    | `make_column_excitation_file` | ✓ | ✓ |
    | `make_Ie_top_file` | ✓ | ✓ |
    | `make_current_file` | ✓ | ✓ |
    | `make_heating_rate_file` | ✓ | ✓ |
    | `make_psd_file` | | ✓ |
    | `downsampling_fluxes` | | ✓ |

## Excitation rates

```@docs; canonical=false
make_volume_excitation_file
make_column_excitation_file
```

## Boundary condition

```@docs; canonical=false
make_Ie_top_file
```

## Current and heating

```@docs; canonical=false
make_current_file
make_heating_rate_file
```

## Phase space density

```@docs; canonical=false
make_psd_file
```

## Downsampling

```@docs; canonical=false
downsampling_fluxes
```
