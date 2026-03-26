# Input Flux

The input flux describes the precipitating electron flux injected at the top of the
simulated ionosphere. It is the primary "forcing" of an AURORA simulation and is fully
specified by an [`InputFlux`](@ref) object, which combines:

- an **energy spectrum** ([`AbstractSpectrum`](@ref)) — the shape of the differential
  number flux as a function of energy, and
- a **temporal modulation** ([`AbstractModulation`](@ref)) — how that flux varies in time.

Once an `InputFlux` is built, it is bundled with an [`AuroraModel`](@ref) into an
[`AuroraSimulation`](@ref) and executed with [`run!`](@ref). Internally, `run!` calls
[`compute_flux`](@ref) to produce the full 3-D flux array `Ie_top` of shape
`[n_beams, n_t, n_E]` (units: #e⁻/m²/s).

!!! note "Units"
    Energy flux (`IeE_tot`) is always given in **W/m²**. Energies and characteristic
    energies (`E₀`, `ΔE`, `E_min`) are in **eV**.

---

## Spectrum types

A spectrum controls the *energy shape* of the input flux. All built-in spectrum types
are sub-types of [`AbstractSpectrum`](@ref) and are normalized so that the integrated
energy flux equals `IeE_tot` (W/m²).

| Type | Shape | Key parameters |
|------|-------|----------------|
| [`FlatSpectrum`](@ref) | Uniform Φ(E) = const above `E_min` | `IeE_tot`, `E_min` |
| [`GaussianSpectrum`](@ref) | Φ(E) ∝ exp(-(E-E₀)²/ΔE²) | `IeE_tot`, `E₀`, `ΔE` |
| [`MaxwellianSpectrum`](@ref) | Maxwellian + optional low-energy tail | `IeE_tot`, `E₀` |
| [`FileSpectrum`](@ref) | Loaded from a `.mat` file | `filename`, `interpolation` |

`FileSpectrum` is the exception: it loads the *complete* flux distribution (including any
time dependence) directly from a file, so it must be used with
[`ConstantModulation`](@ref) and cannot be combined with other modulation types.

```@docs; canonical=false
AbstractSpectrum
FlatSpectrum
GaussianSpectrum
MaxwellianSpectrum
FileSpectrum
evaluate_spectrum
```

---

## Modulation types

A modulation controls the *time dependence* of the flux. All modulation types are
sub-types of [`AbstractModulation`](@ref). The modulation factor (a value between 0 and 1)
is multiplied by the base spectrum at every time step.

| Type | Behaviour |
|------|-----------|
| [`ConstantModulation`](@ref) | Flux is constant (switches on at t = 0). Default. |
| [`SinusoidalFlickering`](@ref) | Smooth oscillation at frequency `f` Hz |
| [`SquareFlickering`](@ref) | On/off square-wave at frequency `f` Hz |
| [`SmoothOnset`](@ref) | C∞-smooth ramp from 0 to 1 over `[t_start, t_end]` |

```@docs; canonical=false
AbstractModulation
ConstantModulation
SinusoidalFlickering
SquareFlickering
SmoothOnset
apply_modulation
```

---

## Building an InputFlux

[`InputFlux`](@ref) is the top-level object that bundles a spectrum, a modulation, the
set of pitch-angle beams that carry the precipitating electrons (`beams`), and an optional
source altitude (`z_source`).

**`beams`** selects which pitch-angle beams receive precipitating flux. Beam indices
correspond to the downward-pointing elements of the pitch-angle grid (i.e. beams with
μ < 0). Upward-going beams are always forced to zero inside `compute_flux`.
The total spectrum is divided among the selected beams weighted by their solid angle
(Ω_beam), so the total energy flux is preserved regardless of how many beams are chosen.

**`z_source`** (km) is the altitude from which electrons are launched. If set, electrons
of different energies and pitch angles acquire a *travel-time delay* before reaching the
top of the ionosphere — faster (higher energy) electrons arrive earlier. This produces
realistic energy- and angle-dependent timing when combined with a flickering modulation.
If omitted (`NaN`, the default), all electrons arrive at the top of the ionosphere
simultaneously with no delay.

### Steady-state example

```julia
# Flat spectrum from 2900 eV to E_max.
# No modulation is needed; ConstantModulation is the default.
flux = InputFlux(FlatSpectrum(1e-2; E_min=2900); beams=1:2)
```

### Time-dependent examples

```julia
# Sinusoidal flickering at 5 Hz, electrons sourced 3000 km above the ionosphere top.
# The travel-time delay is applied automatically for each energy and pitch angle.
flux = InputFlux(FlatSpectrum(1e-2; E_min=100), SinusoidalFlickering(5.0);
                 beams=1, z_source=3000.0)

# Maxwellian spectrum with a smooth onset over the first 0.1 s.
flux = InputFlux(MaxwellianSpectrum(1e-2, 1000.0), SmoothOnset(0.0, 0.1);
                 beams=1:2)

# Square-wave flickering at 10 Hz with 50 % modulation depth.
flux = InputFlux(GaussianSpectrum(1e-2, 5000.0, 500.0), SquareFlickering(10.0; amplitude=0.5);
                 beams=1)
```

### File-based example

```julia
# Load a time-dependent flux directly from a .mat file.
# The file already contains the full temporal behaviour, so no modulation is used.
flux = InputFlux(FileSpectrum("path/to/flux.mat"; interpolation=:linear))
```

```@docs; canonical=false
InputFlux
compute_flux
```

---

## Load from file

[`Ie_top_from_file`](@ref) is the lower-level function called internally when a
[`FileSpectrum`](@ref) is used. It handles mapping the file's time grid onto the
simulation time grid via the chosen interpolation scheme.

### File format requirements

The `.mat` file must contain the variable `Ie_total` with shape **`[n_μ, n_t_file, n_E]`**
(units: #e⁻/m²/s), and optionally `t_top` (a vector of length `n_t_file`, in seconds).

| Dimension | Meaning | Constraint |
|-----------|---------|------------|
| dim 1 — `n_μ` | pitch-angle beams | must **exactly** match the model's beam count |
| dim 2 — `n_t_file` | time steps in the file | flexible (see time-grid notes below) |
| dim 3 — `n_E` | energy bins | must be **≥** the model's `n_E`; only the first `n_E` bins are used |

!!! warning "Strict beam/energy alignment"
    AURORA does **not** interpolate in the pitch-angle or energy dimensions when loading
    from a file. The array is indexed directly, so the grids of the file and the simulation
    model must correspond exactly. If they do not, the results will be silently wrong or an
    error will be raised.

### Preparing the file — matching the grids

Before writing a file, build an `AuroraModel` with the same settings you intend to use
for the simulation, then read the grids off it:

```julia
using AURORA, MAT

# Build the model (same settings as the simulation that will use the file)
θ_lims = 180:-10:0
E_max  = 1000.0
model  = AuroraModel([100, 600], θ_lims, E_max, msis_file, iri_file, 13)

# Extract the grids your file must be aligned to
μ_center = model.pitch_angle_grid.μ_center  # length n_μ
E_centers = model.energy_grid.E_centers     # length n_E
n_μ = length(μ_center)
n_E = length(E_centers)
```

Now build your flux array following the `[n_μ, n_t_file, n_E]` layout.
Downward-going beams correspond to **negative** `μ_center` values; upward-going beams
(positive `μ_center`) will be forced to zero by `Ie_top_from_file` regardless of what
the file contains.

With the standard `θ_lims = 180:-10:0` convention, `μ_center` starts at −1 (straight
down, θ = 180°) and increases toward +1 (straight up, θ = 0°), so **downward-going
beams are the first indices** of the array.

### Writing the file

```julia
# Example: time-dependent flux, 100 time steps
t_file = collect(0.0:0.01:0.99)   # 100 steps, dt = 10 ms
n_t_file = length(t_file)

Ie_total = zeros(n_μ, n_t_file, n_E)
# ... fill Ie_total with your flux values ...

matwrite("my_flux.mat", Dict(
    "Ie_total" => Ie_total,
    "t_top"    => t_file,
))
```

For a steady-state (no time dependence), simply use a single time slice and omit `t_top`:

```julia
Ie_total = zeros(n_μ, 1, n_E)
# ... fill steady-state values ...
matwrite("my_flux_steady.mat", Dict("Ie_total" => Ie_total))
```

### Time-grid handling

If `t_top` is present, the file's time axis does **not** need to match the simulation
time grid in length or step size. The file is always re-sampled:

- The file time is treated as **relative** — `t_top[1]` is always aligned with the
  simulation's `t[1]`, regardless of what the absolute values in `t_top` are. Only
  durations and spacings matter.
- If the file is **shorter** than the simulation, the flux is set to zero for any
  simulation time beyond the file's duration (a warning is emitted).
- If the file is **finer** than the simulation (smaller `dt`), the flux is
  point-sampled, which may miss rapid variations — choose a denser simulation grid
  or a smoother file if this matters.
- The `:constant`, `:linear`, and `:pchip` interpolation options control how the
  file values are mapped onto the simulation grid.

### Using the file in a simulation

```julia
flux = InputFlux(FileSpectrum("my_flux.mat"))               # constant (step) interpolation
flux = InputFlux(FileSpectrum("my_flux.mat"; interpolation=:linear))   # linear interpolation
```

No modulation type should be specified — `ConstantModulation` is both the default and
the only permitted choice for `FileSpectrum`.

```@docs; canonical=false
Ie_top_from_file
```
