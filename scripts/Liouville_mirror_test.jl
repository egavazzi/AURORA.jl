# Collisionless Liouville test for the magnetic mirror force.
#
# A monodirectional beam around 110°-105° pitch angle is injected at the top of a
# collisionless domain with a dipole field, on four pitch-angle grids: uniform in θ and
# uniform in μ = cos(θ), at two beam counts (36 and 72). This compares the same beam budget
# distributed differently: uniform-μ beams have equal solid angles and are automatically
# finest in θ around 90°. For each case the script checks against the adiabatic theory
# (conservation of E(1-μ²)/B with B ∝ (R+h)⁻³):
#   1. drain rate: exponential decay of the uppermost input beam below the top vs the
#      discrete operator prediction  d(ln Ie)/dz = -κ(c_loss + μ̄)/μ̄
#   2. mirror altitude: half-decay altitude of the downgoing flux vs the adiabatic
#      mirror window and its flux-weighted median
#   3. reflection and bottom penetration, and their convergence under μ-refinement
#
# Run with:  julia --project=. scripts/Liouville_mirror_test.jl [dz_max]
# The optional dz_max argument (km, default 25) refines the altitude grid; the drain-rate
# truncation error shrinks linearly with the grid spacing near the top.
using AURORA
using Printf

const R_earth = 6.371e6
dz_max = isempty(ARGS) ? 25.0 : parse(Float64, ARGS[1])

function run_case(θ_lims, input_beams, label; dz_max = 25.0)
    msis_file = find_msis_file()
    iri_file  = find_iri_file()
    model = AuroraModel([100, 600], θ_lims, 100, msis_file, iri_file;
                        magnetic_field = DipoleField(; R_earth))
    model.altitude_grid = AltitudeGrid(100, 600; dz_max)
    flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = input_beams)
    dir = mktempdir()
    sim = AuroraSimulation(model, flux, dir; mode = SteadyStateMode())
    initialize!(sim; verbose = false)

    # Make the run collisionless: no neutrals (elastic + inelastic collisions), no thermal
    # electrons (e-e losses), and no artificial time-of-flight diffusion
    for sp in model.species
        sp.density .= 0.0
    end
    model.ionosphere.ne .= 0.0
    sim.cache.matrices.D .= 0.0

    run!(sim)
    res = load_results(dir)

    # Sum over energy bins (collisionless: all bins behave identically)
    Ie = dropdims(sum(res.Ie[:, :, 1, :]; dims = 3); dims = 3)   # [n_z, n_beams]
    h  = res.h_atm
    n_z = length(h)
    grid = model.pitch_angle_grid
    μ̄ = grid.μ_center
    μ_lims = grid.μ_lims
    down = findall(<(0), μ̄)
    up   = findall(>(0), μ̄)

    down_flux(i) = sum(abs.(μ̄[down]) .* Ie[i, down])
    up_flux(i)   = sum(μ̄[up] .* Ie[i, up])

    # Adiabatic mirror altitudes for the injected pitch-angle range
    h_top = h[end]
    mirror_h(absμ) = (R_earth + h_top) * (1 - absμ^2)^(1/3) - R_earth
    μ_a = abs(μ_lims[last(input_beams) + 1])   # shallow edge -> mirrors high
    μ_b = abs(μ_lims[first(input_beams)])      # steep edge -> mirrors deep
    μ_med = sqrt((μ_a^2 + μ_b^2) / 2)          # flux-weighted median pitch angle

    # Half-decay altitude of the total downgoing flux
    P = [down_flux(i) for i in 1:n_z] ./ down_flux(n_z)
    i2 = findfirst(>(0.5), P)
    f = (0.5 - P[i2 - 1]) / (P[i2] - P[i2 - 1])
    h_half = h[i2 - 1] + f * (h[i2] - h[i2 - 1])

    # Drain rate of the uppermost input beam (its upwind neighbour is empty -> pure decay)
    b0 = first(input_beams)
    Ω = beam_weight(grid.θ_lims)
    c_loss = π * (1 - μ_lims[b0 + 1]^2) / Ω[b0]
    fit_lo, fit_hi = h_top - 100e3, h_top - 5e3
    sel = [i for i in 1:n_z if fit_lo <= h[i] <= fit_hi && Ie[i, b0] > 0]
    x = h[sel]; y = log.(Ie[sel, b0])
    slope_fit = sum((x .- sum(x)/length(x)) .* (y .- sum(y)/length(y))) /
                sum((x .- sum(x)/length(x)).^2)
    κ_mid = 3 / (R_earth + (fit_lo + fit_hi)/2)
    slope_pred = -κ_mid * (c_loss + μ̄[b0]) / μ̄[b0]

    @printf("--- %s (%d beams, input θ = %.1f°-%.1f°) ---\n", label, grid.n_beams,
            acosd(μ_lims[first(input_beams)]), acosd(μ_lims[last(input_beams) + 1]))
    @printf("Drain rate:       %.4e /m simulated vs %.4e /m predicted (%+.1f%%)\n",
            slope_fit, slope_pred, 100 * (slope_fit/slope_pred - 1))
    @printf("Mirror window:    %.1f - %.1f km (adiabatic), median %.1f km\n",
            mirror_h(μ_b)/1e3, mirror_h(μ_a)/1e3, mirror_h(μ_med)/1e3)
    @printf("Half-decay:       %.1f km simulated (%+.1f km from analytic median)\n",
            h_half/1e3, (h_half - mirror_h(μ_med))/1e3)
    @printf("Reflection:       up/down at top = %.4f\n", up_flux(n_z)/down_flux(n_z))
    @printf("Penetration:      down at bottom / top = %.2e\n\n", down_flux(2)/down_flux(n_z))

    return (; slope_err = slope_fit/slope_pred - 1,
              h_half_err = h_half - mirror_h(μ_med),
              window = (mirror_h(μ_b), mirror_h(μ_a)),
              h_half, penetration = down_flux(2)/down_flux(n_z))
end

@printf("=== Collisionless Liouville mirror test (dz_max = %.1f km) ===\n\n", dz_max)

## 36 beams — uniform in θ (5° everywhere) vs uniform in μ (Δμ = 1/18)
θ36 = run_case(180:-5:0, 15:15, "uniform in θ: 5° beams"; dz_max)
μ36 = run_case(acosd.(range(-1, 1, length = 37)), 13:13, "uniform in μ: Δμ = 1/18"; dz_max)

## 72 beams — uniform in θ (2.5° everywhere) vs uniform in μ (Δμ = 1/36)
θ72 = run_case(180:-2.5:0, 29:30, "uniform in θ: 2.5° beams"; dz_max)
μ72 = run_case(acosd.(range(-1, 1, length = 73)), 25:26, "uniform in μ: Δμ = 1/36"; dz_max)

## Verdict
# The mirror term is an advection in μ, so the discretization error is first order in the
# beam width Δμ along the cascade path (110° -> 90°), not in the θ-width near 90°.
# Uniform-θ grids waste μ-resolution near the field-aligned directions (Δμ = sinθ·Δθ), so
# at equal beam count the uniform-μ grid is finer on the cascade:
# θ36 Δμ ≈ 0.085, μ36 Δμ ≈ 0.056, θ72 Δμ ≈ 0.043, μ72 Δμ ≈ 0.028.
# Convergence is therefore expected in that order.
# The drain-rate residual is first-order spatial truncation of the upwind streaming: it
# measures ≈ (drain rate)·Δz/2 across all cases, growing as finer beams drain over fewer
# grid points (e-folding 131 km -> 46 km against Δz ≈ 10 km near the top on the default
# grid). The 15% bound accommodates that truncation at Δμ = 1/36 with dz_max = 25; it
# shrinks linearly when the script is run with a smaller dz_max.
cases = (θ36, μ36, θ72, μ72)
println("=== Verdict ===")
@printf("Operator drain rate within theory:  %s (expect < 15%%)\n",
        join((@sprintf("%.1f%%", 100abs(c.slope_err)) for c in cases), " / "))
@printf("Half-decay inside adiabatic window: %s\n",
        join((string(c.window[1] < c.h_half < c.window[2]) for c in cases), " / "))
@printf("Convergence in Δμ (0.085 -> 0.056 -> 0.043 -> 0.028):\n")
@printf("  |Δh_half|:   %s km\n",
        join((@sprintf("%.1f", abs(c.h_half_err)/1e3) for c in cases), " -> "))
@printf("  penetration: %s\n",
        join((@sprintf("%.3f", c.penetration) for c in cases), " -> "))

ok = all(abs(c.slope_err) < 0.15 for c in cases) &&
     all(c.window[1] < c.h_half < c.window[2] for c in cases) &&
     issorted([abs(c.h_half_err) for c in cases]; rev = true) &&
     issorted([c.penetration for c in cases]; rev = true)
println(ok ? "\nLIOUVILLE TEST PASSED" : "\nLIOUVILLE TEST FAILED")
