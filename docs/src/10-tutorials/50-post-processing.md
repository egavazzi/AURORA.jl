# [Post-Processing & Analysis](@id Post-Processing)

After running a simulation, AURORA provides several post-processing functions that compute
derived physical quantities from the raw electron flux output. Each reads `simulation_data.nc`
and writes its result as a NetCDF file into the `analysis/` subdirectory of the simulation's
save directory — see [Output & data](@ref Output) for the file layout and schema.

This page provides a small description of each analysis function.

!!! tip
    All analysis functions accept either a directory path or an
    [`AuroraSimulation`](@ref) object directly. See the [API reference](@ref Analysis) for
    full docstrings.

## Excitation rates

### Volume excitation rates

[`make_volume_excitation_file`](@ref) reads the electron flux from `simulation_data.nc`,
sums over pitch-angle beams to get the omnidirectional flux, and computes the volume
excitation rate $Q = I_e \times \sigma \times n$ for several optical emissions (4278 Å,
6730 Å, 7774 Å, 8446 Å, O¹D, O¹S) and ionizations (O⁺, O₂⁺, N₂⁺). Results are saved as
`analysis/volume_excitation.nc`.

```julia
make_volume_excitation_file(sim)
```

### Column-integrated excitation rates

[`make_column_excitation_file`](@ref) integrates the volume excitation rates in altitude,
taking into account the finite speed of light (photon travel time from emission altitude to
the bottom of the column). This must be run **after** `make_volume_excitation_file`.
Results are saved as `analysis/column_excitation.nc`.

```julia
make_column_excitation_file(sim)
```

## Top-of-model flux

[`make_Ie_top_file`](@ref) extracts the electron flux at the top of the ionosphere
(maximum simulation altitude), for each pitch-angle beam. This contains *all* beams,
including the upward (backscattered) flux, so it differs from the `Ie_input` boundary
condition stored in `simulation_data.nc`. It is useful for comparing with satellite or rocket
observations. Results are saved as `analysis/Ie_top.nc` (load with [`load_Ie_top`](@ref)).

```julia
make_Ie_top_file(sim)
```

## Field-aligned currents

[`make_current_file`](@ref) computes the field-aligned electron current density $J$ (A/m²) and
electron energy flux (eV/m²/s), separated into upward and downward components. Results are saved
as `analysis/currents.nc`.

```julia
make_current_file(sim)
```

## Heating rates

[`make_heating_rate_file`](@ref) computes the rate at which energy is transferred from
superthermal (energetic) electrons to thermal electrons through Coulomb collisions.
Results are saved as `analysis/heating_rate.nc`.

```julia
make_heating_rate_file(sim)
```

## Phase-space density

[`make_psd_file`](@ref) converts the electron flux into phase-space density. It can
compute:
- the full distribution `f(z, μ, t, E)` (`:f_only`),
- the reduced parallel distribution `F(z, t, v∥)` (`:F_only`), or
- both (`:both`).

Results are saved as `analysis/psd.nc`.

```julia
make_psd_file(sim; compute = :both)
```

You can also pass custom `v_parallel` bin edges:

```julia
make_psd_file(sim; compute = :F_only, vpar_edges = range(-2e7, 2e7; length = 100))
```