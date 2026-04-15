# [Post-Processing & Analysis](@id Post-Processing)

After running a simulation, AURORA provides several post-processing functions that compute
derived physical quantities from the raw electron flux output and save them alongside the
simulation data.

This tutorial walks through each analysis function in detail. We first run a short
simulation to produce output files, then demonstrate each post-processing step.

!!! tip
    All analysis functions accept either a directory path or an
    [`AuroraSimulation`](@ref) object directly. See the [API reference](@ref Analysis) for
    full docstrings.

## Setup: run a short simulation

```@example post_proc
using AURORA

msis_file = find_msis_file()
iri_file  = find_iri_file()

model = AuroraModel(
    [100, 600],    # altitude limits [km]
    180:-45:0,     # pitch-angle bin edges [°]
    1000,          # maximum energy [eV]
    msis_file,
    iri_file,
    13             # magnetic field angle to zenith [°]
)

flux = InputFlux(
    FlatSpectrum(1e-3; E_min=100),
    SinusoidalFlickering(5.0);
    beams=1,
    z_source=3000.0,
)

savedir = mkpath(joinpath("data", "post_proc_tutorial"))
savedir = mktempdir()  # hide — redirect to OS temp so .mat files are not deployed to gh-pages

sim = AuroraSimulation(
    model, flux,
    0.2, 0.01, savedir;
    CFL_number=128,
    max_memory_gb=4.0,
)

run!(sim)
nothing # hide
```

## Excitation rates

### Volume excitation rates

[`make_volume_excitation_file`](@ref) loads the electron flux from each
`IeFlickering-NN.mat` file, sums over pitch-angle beams to get the omnidirectional flux,
and computes the volume excitation rate $Q = I_e \times \sigma \times n$ for several
optical emissions (4278 Å, 6730 Å, 7774 Å, 8446 Å, O¹D, O¹S) and ionizations (O⁺, O₂⁺,
N₂⁺). Results are saved as `Qzt_all_L.mat`.

```@example post_proc
make_volume_excitation_file(sim)
nothing # hide
```

### Column-integrated excitation rates

[`make_column_excitation_file`](@ref) integrates the volume excitation rates in altitude,
taking into account the finite speed of light (photon travel time from emission altitude to
the bottom of the column). This must be run **after** `make_volume_excitation_file`.
Results are saved as `I_lambda_of_t.mat`.

```@example post_proc
make_column_excitation_file(sim)
nothing # hide
```

## Boundary condition (top flux)

[`make_Ie_top_file`](@ref) extracts the electron flux at the top of the ionosphere
(maximum simulation altitude), for each pitch-angle beam. This is useful for verifying the
applied boundary condition and for comparing with satellite or rocket observations.
Results are saved as `Ie_top.mat`.

```@example post_proc
make_Ie_top_file(sim)
nothing # hide
```

## Field-aligned currents

[`make_current_file`](@ref) computes the field-aligned electron current density $J$ (A/m²) and
electron energy flux (eV/m²/s), separated into upward and downward components. Results are saved
as `J.mat`.

```@example post_proc
make_current_file(sim)
nothing # hide
```

## Heating rates

[`make_heating_rate_file`](@ref) computes the rate at which energy is transferred from
superthermal (energetic) electrons to thermal electrons through Coulomb collisions.
Results are saved as `heating_rate.mat`.

```@example post_proc
make_heating_rate_file(sim)
nothing # hide
```

## Phase-space density

[`make_psd_file`](@ref) converts the electron flux into phase-space density. It can
compute:
- the full distribution `f(z, μ, t, E)` (`:f_only`),
- the reduced parallel distribution `F(z, t, v∥)` (`:F_only`), or
- both (`:both`).

```@example post_proc
make_psd_file(sim; compute = :both)
nothing # hide
```

You can also pass custom `v_parallel` bin edges:

```julia
make_psd_file(sim; compute = :F_only, vpar_edges = range(-2e7, 2e7; length = 100))
```

## Downsampling

[`downsampling_fluxes`](@ref) reduces the time resolution of the saved electron flux files.
This is useful when transferring data between machines or when finer time resolution is not
needed for a given analysis. The downsampled files are saved in a subdirectory named
`downsampled_Xx` (where `X` is the downsampling factor).

```@example post_proc
downsampling_fluxes(sim, 2)
nothing # hide
```

## Output summary

After running all post-processing steps, the output directory contains:

```@example post_proc
readdir(savedir)
```

See the [Analysis API reference](@ref Analysis) for complete function documentation.
