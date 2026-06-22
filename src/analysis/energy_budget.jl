# ======================================================================================== #
#                              ENERGY-BUDGET DIAGNOSTIC                                     #
# ======================================================================================== #
#
# Global steady-state energy balance for the suprathermal electron population, in the spirit
# of the conservation check that the TRANS family (transcar / aeroplanets) prints. At steady
# state, the downward energy flux entering the top must be accounted for by:
#
#     E_in  =  E_inelastic  +  E_heating  +  E_escape  +  residual
#
# where
#   E_in        vertical energy flux carried by the *downward* beams at the top boundary
#               (= the imposed precipitation), Σ_beam Σ_E Ie·E·|μ|   [eV m⁻² s⁻¹]
#   E_escape    vertical energy flux carried by the *upward* beams at the top boundary
#               (backscattered electrons leaving the domain)
#   E_inelastic energy deposited in neutral excitation + ionization thresholds,
#               Σ_z Δz Σ_species n(z) Σ_levels threshold·Σ_E Ie_omni·σ
#   E_heating   energy transferred to the thermal electrons (Coulomb), via
#               `calculate_heating_rate`, integrated over altitude
#   residual    everything unaccounted: physically the energy of electrons degraded below the
#               grid floor (sub-floor thermalisation, a small *positive* term on a good grid),
#               PLUS any numerical non-conservation. A healthy run shows a small positive
#               residual; a grid that violates the ΔE < threshold stability bound shows the
#               inelastic term ballooning and the residual going large/negative (energy
#               "created") — exactly the failure this diagnostic is meant to catch.
#
# The vertical-energy-flux projection (Σ Ie·E·|μ|) matches `field_aligned_beam_norm`, so
# E_in computed here equals the InputFlux's `IeE_tot` normalisation.

"""
    EnergyBudget

Result of [`energy_budget`](@ref). All energy fluxes are column-integrated and in eV m⁻² s⁻¹.

Fields: `input`, `inelastic`, `heating`, `escape`, `residual`, and `residual_fraction`
(= `residual / input`). `inelastic_by_species` breaks the inelastic term down per species.
"""
struct EnergyBudget
    input::Float64
    inelastic::Float64
    heating::Float64
    escape::Float64
    residual::Float64
    residual_fraction::Float64
    inelastic_by_species::Vector{Pair{String,Float64}}
end

function Base.show(io::IO, ::MIME"text/plain", b::EnergyBudget)
    pct(x) = round(100x / b.input, digits=2)
    println(io, "EnergyBudget (eV m⁻² s⁻¹, steady state):")
    println(io, "├── input (precip., ↓ top) : ", b.input)
    println(io, "├── inelastic deposited     : ", b.inelastic, "  (", pct(b.inelastic), "%)")
    for (name, val) in b.inelastic_by_species
        println(io, "│     ├── ", rpad(name, 4), "             : ", val, "  (", pct(val), "%)")
    end
    println(io, "├── thermal heating         : ", b.heating, "  (", pct(b.heating), "%)")
    println(io, "├── escape (↑ top)          : ", b.escape, "  (", pct(b.escape), "%)")
    println(io, "└── residual (subfloor+num) : ", b.residual, "  (", b.residual_fraction * 100, "%)")
end

"""
    energy_budget(sim::AuroraSimulation; verbose=true) -> EnergyBudget

Compute the steady-state energy balance of a finished simulation directly from the in-memory
solution (`sim.cache.Ie`). Returns an [`EnergyBudget`](@ref); also prints a summary unless
`verbose=false`.

Use as a guardrail: on a converged, stable run the residual is a small positive fraction (the
sub-floor thermalisation that AURORA does not track on-grid). A large or negative residual
flags energy non-conservation — e.g. a grid whose maximum bin width exceeds the lowest
ionization threshold, which destabilises the high→low energy-degradation sweep.
"""
function energy_budget(sim::AuroraSimulation; verbose::Bool = true)
    model = sim.model
    eg    = model.energy_grid
    z     = model.altitude_grid.h          # m
    μ     = model.pitch_angle_grid.μ_center
    E     = eg.E_centers                   # eV
    ne    = model.ionosphere.ne
    Te    = model.ionosphere.Te
    n_z   = length(z)
    n_μ   = length(μ)
    n_E   = length(E)

    # Steady-state snapshot of the flux: reshape [n_z·n_μ, n_t, n_E] → [n_z, n_μ, n_E]
    # (row = (i_μ-1)·n_z + i_z). Use the last time index (n_t = 1 for steady state).
    Ie_raw = sim.cache.Ie
    it = size(Ie_raw, 2)
    Ie = reshape(@view(Ie_raw[:, it, :]), n_z, n_μ, n_E)   # [n_z, n_μ, n_E]

    # ---- Boundary energy fluxes at the top altitude (i_z = n_z) --------------------------
    # Vertical energy flux = Σ_beam Σ_E Ie·E·|μ| (matches field_aligned_beam_norm).
    input  = 0.0   # downward beams (μ < 0): the imposed precipitation
    escape = 0.0   # upward beams (μ > 0): backscattered electrons leaving the top
    @inbounds for iμ in 1:n_μ, iE in 1:n_E
        f = Ie[n_z, iμ, iE] * E[iE] * abs(μ[iμ])
        μ[iμ] < 0 ? (input += f) : (escape += f)
    end

    # ---- Omnidirectional flux for the volumetric (deposition / heating) terms -----------
    Ie_omni = dropdims(sum(Ie, dims=2), dims=2)            # [n_z, n_E]

    # ---- Inelastic energy deposition: Σ_z Δz Σ_sp n Σ_levels threshold·Σ_E Ie_omni·σ -----
    dep_profile_total = zeros(n_z)
    inelastic_by_species = Pair{String,Float64}[]
    for sp in model.species
        σ      = sp.cross_sections          # [n_levels, n_E]
        levels = sp.excitation_levels       # [n_levels, 2]: (energy loss, #secondaries)
        dens   = sp.density                 # [n_z]
        dep_sp = zeros(n_z)
        for lvl in axes(levels, 1)
            E_loss = levels[lvl, 1]
            E_loss > 0 || continue          # skip elastic (no energy loss)
            for iz in 1:n_z
                acc = 0.0
                @inbounds for iE in 1:n_E
                    acc += Ie_omni[iz, iE] * σ[lvl, iE]
                end
                dep_sp[iz] += dens[iz] * acc * E_loss
            end
        end
        col = trapz_altitude(z, dep_sp)
        push!(inelastic_by_species, String(sp.name) => col)
        dep_profile_total .+= dep_sp
    end
    inelastic = trapz_altitude(z, dep_profile_total)

    # ---- Thermal-electron heating (reuse the existing Coulomb-loss routine) --------------
    heating_profile = calculate_heating_rate(z, [0.0], reshape(Ie_omni, n_z, 1, n_E),
                                             E, ne, Te)[:, 1]
    heating = trapz_altitude(z, heating_profile)

    # ---- Closure ------------------------------------------------------------------------
    residual = input - inelastic - heating - escape
    budget = EnergyBudget(input, inelastic, heating, escape, residual,
                          input > 0 ? residual / input : NaN, inelastic_by_species)
    verbose && show(stdout, MIME"text/plain"(), budget)
    verbose && println()
    return budget
end

# Simple trapezoidal integral over a (possibly non-monotone / descending) altitude grid.
function trapz_altitude(z, f)
    s = 0.0
    @inbounds for i in 1:(length(z) - 1)
        s += 0.5 * (f[i] + f[i + 1]) * abs(z[i + 1] - z[i])
    end
    return s
end
