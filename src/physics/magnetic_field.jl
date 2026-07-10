"""
    AbstractFieldProfile

Supertype for magnetic field strength profiles used by the mirror force.

A concrete subtype is a callable `field(h)` returning the (relative) magnetic field strength
at altitude `h` in meters. Only the logarithmic gradient of the field enters the mirror force,
so any normalization is fine — `field(h) ∝ B(h)` is enough.

The mirror-force operator needs the inverse focusing length `κ = −d(ln B)/ds` along the field
line. By default this is obtained by finite-differencing `log.(field.(h))` along
`s = h / cosd(B_angle_to_zenith)` (see [`mirror_kappa`](@ref)). A subtype may optionally
specialize [`mirror_kappa`](@ref) to return `κ` in closed form, as [`DipoleField`](@ref) does.

Pass an instance to [`AuroraModel`](@ref) with the `magnetic_field` keyword to enable the
magnetic mirror force.
"""
abstract type AbstractFieldProfile end

"""
    DipoleField(; R_earth = 6.371e6)

Relative magnetic field strength along a dipole field line, as a function of altitude in
meters: `B(h) ∝ (R_earth / (R_earth + h))³` for a radial line. Only the logarithmic gradient
of the field matters for the mirror force, so the profile is normalized to `B(0) = 1`.

Pass it to [`AuroraModel`](@ref) with the `magnetic_field` keyword to enable the magnetic
mirror force:

```julia
model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file;
                    magnetic_field = DipoleField())
```

The inverse focusing length `κ = −d(ln B)/ds` is computed in closed form for the inclined
dipole field line whose dip matches `B_angle_to_zenith` (see [`mirror_kappa`](@ref)),
including the latitude-dependent dipole factor — it is *not* obtained by differencing the
radial `B(h)` profile above. At `B_angle_to_zenith = 0` the line is radial and `κ = 3/r`.

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
struct DipoleField <: AbstractFieldProfile
    R_earth::Float64
end

DipoleField(; R_earth = 6.371e6) = DipoleField(R_earth)

(field::DipoleField)(h) = (field.R_earth / (field.R_earth + h))^3

"""
    mirror_kappa(field, h, B_angle_to_zenith)

Inverse focusing length `κ = −d(ln B)/ds` [1/m] at every altitude in `h` [m], along the field
line inclined by `B_angle_to_zenith` degrees from the zenith (`s = h / cosd(angle)`).

The default method for any [`AbstractFieldProfile`](@ref) finite-differences `log.(field.(h))`
along `s` (one-sided at the domain edges, central in the interior). The upwind mirror-force
discretization assumes a field converging downward, so a negative `κ` (field strengthening
with altitude) is rejected with an `ArgumentError`.

Concrete profiles may specialize this method to return `κ` in closed form; see the
[`DipoleField`](@ref) method below.
"""
function mirror_kappa(field::AbstractFieldProfile, h, B_angle_to_zenith)
    s = h ./ cosd(B_angle_to_zenith)
    lnB = log.(field.(h))
    n_z = length(s)
    κ = zeros(n_z)
    κ[1] = -(lnB[2] - lnB[1]) / (s[2] - s[1])
    for iz in 2:(n_z - 1)
        κ[iz] = -(lnB[iz + 1] - lnB[iz - 1]) / (s[iz + 1] - s[iz - 1])
    end
    κ[n_z] = -(lnB[n_z] - lnB[n_z - 1]) / (s[n_z] - s[n_z - 1])
    any(<(0), κ) && throw(ArgumentError(
        "the magnetic field strength must be non-increasing with altitude: the upwind \
         mirror-force discretization assumes a field converging downward (μ-drift toward \
         μ = +1 everywhere)."))
    return κ
end

"""
    mirror_kappa(field::DipoleField, h, B_angle_to_zenith)

Closed-form inverse focusing length `κ` [1/m] along the inclined dipole field line whose dip
matches `B_angle_to_zenith`.

For a dipole, `B ∝ √(1 + 3sin²λ) / r³` with `cos²λ = r / r_eq` along the field line, where `λ`
is the magnetic latitude and `r_eq` the equatorial-crossing distance. Converting `d/dr` to
`d/ds` with `ds = dr / cos(angle)` gives

    κ = cos(angle) · [3/r + (3/2) / (r_eq · (1 + 3sin²λ))].

The field line is anchored at the bottom of the domain: the dip relation `tan(I) = 2 tan(λ)`
(with dip `I = 90° − angle`) fixes the reference latitude, and hence `r_eq`, from the bottom
altitude. At `B_angle_to_zenith = 0` the line is radial (`r_eq = ∞`) and `κ = 3/r` exactly.
"""
function mirror_kappa(field::DipoleField, h, B_angle_to_zenith)
    r = field.R_earth .+ h
    cosang = cosd(B_angle_to_zenith)

    # Magnetic latitude of the field line from the dipole dip relation tan(I) = 2 tan(λ),
    # with the dip angle I = 90° − angle measured from the horizontal. Anchored at the
    # bottom of the domain to fix the field line's equatorial-crossing distance r_eq.
    dip = 90.0 - B_angle_to_zenith
    λ_ref = atand(tand(dip) / 2)
    r_eq = r[1] / cosd(λ_ref)^2          # → ∞ (in practice huge) at the pole, λ_ref → 90°

    r_eq > maximum(r) || throw(ArgumentError(
        "The dipole field line implied by B_angle_to_zenith = $(B_angle_to_zenith)° " *
        "crosses the magnetic equator within the altitude domain. " *
        "Use a smaller B_angle_to_zenith or a lower top altitude."))

    κ = similar(r, Float64)
    @inbounds for i in eachindex(r)
        sin2λ = 1.0 - r[i] / r_eq
        κ[i] = cosang * (3.0 / r[i] + 1.5 / (r_eq * (1.0 + 3.0 * sin2λ)))
    end
    return κ
end
