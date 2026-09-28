# ======================================================================================== #
#                 ENERGY OF THE ELECTRONS DEGRADED BELOW THE ENERGY GRID                    #
# ======================================================================================== #
#
# The energy grid starts at a lowest edge E_edges[1] of a few eV. An inelastic collision or
# the Coulomb loss can take an electron below that edge; it then leaves the suprathermal
# population and thermalises. The factors below give the energy those electrons carry, per
# collision (per unit Le·Ie for the Coulomb loss), from the kinematics of each process: the
# arrival range of a non-ionizing channel, and the secondary law of an ionizing channel. They
# do not depend on where the solver places the electrons that stay on the grid, so any
# difference between the energy the solver removes and the energy these factors, the
# inelastic losses and the heating account for is left to the budget's residual.

"""
    subfloor_energy_factors(sp, energy_grid) -> Matrix{Float64}

Energy (eV) carried below the lowest grid edge `E_edges[1]` by the electrons that come out
of one collision of species `sp`. The matrix has the shape of `sp.cross_sections`, row for
row (row 1, elastic, is zero), and one column per energy bin. Multiplied by the collision
rate `n·σ·Ie` and summed, it gives the collisional part of the `subfloor` term of
[`EnergyBudget`](@ref).

For a primary in bin `iE`, of width `ΔE[iE]`, primary energies are taken uniform over the
bin:
- non-ionizing channel of energy loss `L`: the electrons that leave the bin arrive uniformly
  over `[E_edges[iE] - L, min(E_edges[iE+1] - L, E_edges[iE])]`, with density `1/ΔE[iE]`
  per collision. The factor is the energy of the part of that range below `E_edges[1]`;
  arrivals below zero energy carry none.
- ionizing channel of threshold `I`: the energy of the secondaries and of the degraded
  primary that end below `E_edges[1]`, integrated over the secondary law and over the
  primary energies of the bin above `I`, divided by the number of ionization events of the
  bin (`event_count`), the normalization of the cascading matrices. Double ionization uses
  the joint densities of the cascading calculation.
"""
function subfloor_energy_factors(sp, energy_grid)
    E_edges = energy_grid.E_edges
    ΔE = energy_grid.ΔE
    σ = sp.cross_sections
    levels = sp.excitation_levels
    Base.require_one_based_indexing(E_edges, ΔE, σ, levels)
    size(σ) == (size(levels, 1), length(ΔE)) ||
        throw(DimensionMismatch("$(sp.name): cross_sections is $(size(σ)); expected " *
                                "$((size(levels, 1), length(ΔE))) (levels × energy bins)"))

    channel_rows = axes(levels, 1)[2:end]           # row 1 is elastic
    factors = zeros(size(σ))
    for l in channel_rows, iE in eachindex(ΔE)
        levels[l, 2] <= 0 || continue
        factors[l, iE] = excitation_subfloor_energy(E_edges, ΔE, iE, levels[l, 1])
    end

    ionizing_rows = [l for l in channel_rows if levels[l, 2] > 0]
    isempty(ionizing_rows) && return factors
    cache = sp.cascading_data
    cache.E_edges == E_edges ||
        throw(ArgumentError("the cascading matrices of $(sp.name) are on a different " *
                            "energy grid from the model; initialize! the model first"))
    for l in ionizing_rows
        levels[l, 2] in (1, 2) ||
            throw(ArgumentError("$(sp.name): channel row $l ejects $(levels[l, 2]) " *
                                "secondaries; only single and double ionization are modelled"))
    end
    # The law is evaluated millions of times: unwrap it once and enter the integration
    # through a single invokelatest, as calculate_cascading_matrices does.
    law = runtime_law(cache.spec.secondary_law)
    Base.invokelatest(fill_ionization_subfloor_energy!, factors, law, levels, ionizing_rows,
                      E_edges)

    # Per ionization event. A bin whose primary energies all lie below the threshold has no
    # events, and the solver refuses a nonzero cross section there once it ionizes at all.
    min_ionization_E = minimum(levels[l, 1] for l in ionizing_rows)
    for l in ionizing_rows, iE in eachindex(ΔE)
        events = event_count(cache, iE, levels[l, 1])
        if events > 0
            factors[l, iE] /= events
        else
            σ[l, iE] > 0 && min_ionization_E < E_edges[iE] && throw(ArgumentError(
                "$(sp.name): ionizing channel at $(levels[l, 1]) eV has cross section \
                 $(σ[l, iE]) m² in energy bin $(iE) but no ionization events in its \
                 cascading matrices"))
            factors[l, iE] = 0.0
        end
    end
    return factors
end

# Energy below E_edges[1] of the electrons that one collision of energy loss `L` takes out of
# bin `iE`: the leaving range at density 1/ΔE[iE], clipped to [0, E_edges[1]].
function excitation_subfloor_energy(E_edges, ΔE, iE, L)
    lower = max(E_edges[iE] - L, 0.0)
    upper = min(E_edges[iE + 1] - L, E_edges[iE], E_edges[1])
    upper > lower || return 0.0
    return (upper^2 - lower^2) / (2 * ΔE[iE])
end

# Fill `factors[l, iE]`, for the ionizing rows, with the sub-floor energy of the bin's
# ionization events before normalization by the event count.
function fill_ionization_subfloor_energy!(factors, law, levels, ionizing_rows, E_edges)
    tasks = [(l, iE) for l in ionizing_rows for iE in axes(factors, 2)]
    Threads.@threads for k in eachindex(tasks)
        l, iE = tasks[k]
        threshold = levels[l, 1]
        factors[l, iE] = if levels[l, 2] == 1
            single_ionization_subfloor_energy(law, threshold, E_edges[iE], E_edges[iE + 1],
                                              E_edges[1])
        else
            double_ionization_subfloor_energy(law, threshold, E_edges[iE], E_edges[iE + 1],
                                              E_edges)
        end
    end
    return factors
end

# The subintervals of [lower, upper] cut at the `breaks` inside it: the sub-floor integration
# limits have kinks at excess energies that are multiples of the floor energy.
function integration_pieces(lower, upper, breaks)
    points = sort!([lower; [b for b in breaks if lower < b < upper]; upper])
    return [(points[k], points[k + 1]) for k in 1:length(points) - 1]
end

# Single ionization, primary energies E_p ∈ [E_lower, E_upper] above `threshold`, excess
# energy W = E_p - threshold shared as (E_s, W - E_s) with E_s ≤ W/2: the secondary is below
# the floor for E_s < floor, the degraded primary for E_s > W - floor.
function single_ionization_subfloor_energy(law, threshold, E_lower, E_upper, floor_E)
    lower = max(E_lower, threshold)
    total = 0.0
    E_upper > lower || return total
    breaks = (threshold + floor_E, threshold + 2 * floor_E)
    for (a, b) in integration_pieces(lower, E_upper, breaks)
        midpoint, halfwidth = (a + b) / 2, (b - a) / 2
        for (x, w) in zip(GL16_X, GL16_W)
            E_p = midpoint + halfwidth * x
            W = E_p - threshold
            secondary = gauss_legendre16(secondary_energy_density, (law, E_p),
                                         0.0, min(floor_E, W / 2))
            primary = gauss_legendre16(degraded_energy_density, (law, E_p, W),
                                       max(0.0, W - floor_E), W / 2)
            total += halfwidth * w * (secondary + primary)
        end
    end
    return total
end

secondary_energy_density((law, E_p), E_s) = E_s * checked_secondary_law(law, E_s, E_p)
degraded_energy_density((law, E_p, W), E_s) = (W - E_s) * checked_secondary_law(law, E_s, E_p)

# Double ionization: two secondaries with the per-secondary density
# `double_secondary_density` (each at most W/2), and the degraded primary with
# `double_primary_density_cdf` (between W/3 and W).
function double_ionization_subfloor_energy(law, threshold, E_lower, E_upper, E_edges)
    floor_E = E_edges[1]
    lower = max(E_lower, threshold)
    total = 0.0
    E_upper > lower || return total
    breaks = (threshold + floor_E, threshold + 2 * floor_E, threshold + 3 * floor_E)
    for (a, b) in integration_pieces(lower, E_upper, breaks)
        midpoint, halfwidth = (a + b) / 2, (b - a) / 2
        for (x, w) in zip(GL8_X, GL8_W)
            E_p = midpoint + halfwidth * x
            W = E_p - threshold
            W > 0 || continue
            cdf = build_secondary_cdf(E_edges, W / 2, E_p, law)
            secondary = gauss_legendre16(double_secondary_energy_density, (threshold, cdf),
                                         0.0, min(floor_E, W / 2))
            primary = gauss_legendre16(double_primary_energy_density, (threshold, cdf),
                                       W / 3, min(W, floor_E))
            total += halfwidth * w * (2 * secondary + primary)
        end
    end
    return total
end

double_secondary_energy_density(context, E_s) = E_s * double_secondary_density(context, E_s)
double_primary_energy_density(context, E_d) = E_d * double_primary_density_cdf(context, E_d)

"""
    coulomb_subfloor_factors(energy_grid) -> Vector{Float64}

Energy carried below the lowest grid edge by the Coulomb (thermal-electron) loss, per unit
of `Le·Ie`, one value per energy bin. The solver takes `Le·Ie/ΔE[1]` electrons per unit
length out of bin 1 through its lower edge. The loss is continuous, so they cross the floor
at `E_edges[1]` and each carries that energy: the factor is `E_edges[1]/ΔE[1]` for bin 1.
From higher bins the Coulomb loss moves electrons to the bin below, which stays on the grid,
so their factor is zero.
"""
function coulomb_subfloor_factors(energy_grid)
    E_edges = energy_grid.E_edges
    ΔE = energy_grid.ΔE
    Base.require_one_based_indexing(E_edges, ΔE)
    factors = zeros(length(ΔE))
    isempty(factors) || (factors[1] = E_edges[1] / ΔE[1])
    return factors
end

# The per-model inputs of the budget's `subfloor` term: one factor matrix per species, and
# the Coulomb factors times the loss function Le, [n_z, n_E]. They depend only on the model,
# so a time integral computes them once for all its slices.
function budget_subfloor_factors(model)
    eg = model.energy_grid
    collisions = [subfloor_energy_factors(sp, eg) for sp in model.species]
    c = coulomb_subfloor_factors(eg)
    ne = model.ionosphere.ne
    Te = model.ionosphere.Te
    coulomb = zeros(length(ne), length(c))
    for iE in eachindex(c)
        coulomb[:, iE] .= c[iE] .* loss_to_thermal_electrons(eg.E_centers[iE], ne, Te)
    end
    return (; collisions, coulomb)
end
