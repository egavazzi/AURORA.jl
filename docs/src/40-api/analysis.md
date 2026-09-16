# [Analysis](@id Analysis)

Post-processing functions that compute derived quantities from raw simulation output. Each
reads `simulation_data.nc` and writes a file into the `analysis/` subdirectory of the
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
    | `make_energy_budget_file` | ✓ | ✓ |

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

## Energy budget

Global energy balance of the electron population: how much of the precipitating energy flux
ends up in neutral excitation and ionization, in thermal-electron heating, and back out of
the top of the ionosphere.

```@docs; canonical=false
energy_budget
EnergyBudget
make_energy_budget_file
load_energy_budget
```
