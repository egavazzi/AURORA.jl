"""
    DipoleField(; R_earth = 6.371e6)

Relative magnetic field strength along a (near-radial) dipole field line, as a function of
altitude in meters: `B(h) ∝ (R_earth / (R_earth + h))³`. Only the logarithmic gradient of
the field matters for the mirror force, so the profile is normalized to `B(0) = 1`.

Pass it to [`AuroraModel`](@ref) with the `magnetic_field` keyword to enable the magnetic
mirror force:

```julia
model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file;
                    magnetic_field = DipoleField())
```

# Grid resolution for mirroring

The mirror force is an advection in `μ = cos(θ)` (it drifts electrons toward `μ = +1` as the
field converges downward), discretized first order upwind. Two grid choices set its accuracy,
and they must be refined *together*: refining one alone only trades one first-order error for
the other.

- **Pitch angle — use beams uniform in `μ`, not in `θ`.** The truncation error is first order
  in the beam width `Δμ` along the cascade toward `μ = +1`, so `Δμ` is what matters, not `Δθ`.
  A uniform-`θ` grid has `Δμ = sinθ · Δθ`, which is *coarsest at 90°* — exactly where
  near-trapped electrons turn around, wasting resolution on the field-aligned directions.
  Build uniform-`μ` beams straight from `θ_lims`:
  ```julia
  N = 36                                          # number of beams
  θ_lims = acosd.(range(-1, 1, length = N + 1))   # equal Δμ = 2/N, finest near 90°
  ```
  At equal beam count this is strictly better for mirroring: ~30–40 beams put the mirror
  altitudes inside the adiabatic window, and `Δμ ≲ 0.02` (≳ 100 beams) tightens the reflected
  fraction.

- **Altitude — keep `dz_max` small near the top.** Near-90° electrons mirror high in the
  domain, and a beam's drain rate below the injection altitude carries a spatial truncation
  error `≈ (drain rate) · Δz / 2`. Finer beams drain faster, so they demand a finer `Δz` to
  stay resolved. The default spacing law only reaches its `dz_max` cap well above a few
  hundred km, so for a typical 100–600 km domain the top spacing is ~10–12 km *regardless* of
  `dz_max ≥ 10`. Reassign the grid with a smaller cap to actually refine the reflection region:
  ```julia
  model.altitude_grid = AltitudeGrid(100, 600; dz_max = 5)   # km
  ```
"""
struct DipoleField
    R_earth::Float64
end

DipoleField(; R_earth = 6.371e6) = DipoleField(R_earth)

(field::DipoleField)(h) = (field.R_earth / (field.R_earth + h))^3
