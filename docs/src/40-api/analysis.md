# [Analysis](@id Analysis)

Post-processing functions that compute derived quantities from raw simulation output. Each
reads `simulation_data.nc` and writes a NetCDF file into the `analysis/` subdirectory of the
simulation's save directory (see [Output & data](@ref Output) for the layout and schema).

!!! note "Compatibility"

    | Function | Steady-State results | Time-Dependent results |
    |----------|:------------:|:--------------:|
    | `make_volume_excitation_file` | ✓ | ✓ |
    | `make_column_excitation_file` | ✓ | ✓ |
    | `make_Ie_top_file` | ✓ | ✓ |
    | `make_current_file` | ✓ | ✓ |
    | `make_heating_rate_file` | ✓ | ✓ |
    | `make_psd_file` | | ✓ |

## Excitation rates

```@docs; canonical=false
make_volume_excitation_file
make_column_excitation_file
```

## Top-of-model flux

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
