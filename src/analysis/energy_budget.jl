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
#               Σ_s Δs Σ_species n(s) Σ_levels threshold·Σ_E Ie_omni·σ
#   E_heating   energy transferred to the thermal electrons (Coulomb), via
#               `calculate_heating_rate`, integrated along the field line
#   residual    everything unaccounted: physically the energy of electrons degraded below the
#               grid floor (sub-floor thermalisation, a small *positive* term on a good grid),
#               PLUS any numerical non-conservation. A healthy run shows a small positive
#               residual; a grid that violates the ΔE < threshold stability bound shows the
#               inelastic term ballooning and the residual going large/negative (energy
#               "created") — exactly the failure this diagnostic is meant to catch.
#
# Volume integrals run ALONG THE MAGNETIC FIELD LINE (over `model.s_field = h / cos B_angle`),
# NOT over vertical altitude: AURORA's transport conserves the field-aligned flux (μ ∂Ie/∂s),
# so the column path length is s. Integrating over z instead injects a spurious 1/cos(B_angle)
# factor (it vanishes at B_angle = 0, where s == z).
#
# The vertical-energy-flux projection (Σ Ie·E·|μ|) matches `field_aligned_beam_norm`, so
# E_in computed here equals the InputFlux's `IeE_tot` normalisation.

"""
    EnergyBudget

Result of [`energy_budget`](@ref). Energy fluxes are column-integrated along the field line
and in eV m⁻² s⁻¹; reaction rates are in m⁻² s⁻¹; `z_centroid` is in km.

Energy fluxes
- `input`, `escape`   |μ|-weighted *vertical* energy flux of the downward / upward beams at the
  top boundary (`input` matches `IeE_tot`).
- `net`               `input - escape`: the energy actually deposited in the column.
- `inelastic`         energy into neutral excitation + ionization thresholds.
- `ionization`, `excitation`  the inelastic term split into ionizing (≥1 secondary) and
  non-ionizing channels (`inelastic == ionization + excitation`).
- `heating`           energy transferred to the thermal electrons (Coulomb).
- `residual`          `input - inelastic - heating - escape` (sub-floor thermalisation +
  numerical non-conservation); `residual_fraction == residual / input`.
- `albedo`            `escape / input`.
- `input_raw`, `escape_raw`  the *along-field*, un-|μ|-weighted top fluxes, Σ Ie·E.

Rates / geometry
- `ionpairs`, `excevents`   column ion-pair and excitation-event production rates.
- `z_centroid`              energy-deposition-weighted mean altitude (km).
- `inelastic_by_species`    the inelastic term broken down per species.
"""
struct EnergyBudget
    input::Float64
    escape::Float64
    net::Float64
    inelastic::Float64
    ionization::Float64
    excitation::Float64
    heating::Float64
    residual::Float64
    residual_fraction::Float64
    albedo::Float64
    input_raw::Float64
    escape_raw::Float64
    ionpairs::Float64
    excevents::Float64
    z_centroid::Float64
    inelastic_by_species::Vector{Pair{String,Float64}}
end

function Base.show(io::IO, ::MIME"text/plain", b::EnergyBudget)
    pct(x) = round(100x / b.input, digits=2)
    println(io, "EnergyBudget (eV m⁻² s⁻¹, steady state, ∫ along field line):")
    println(io, "├── input (precip., ↓ top) : ", b.input)
    println(io, "├── inelastic deposited     : ", b.inelastic, "  (", pct(b.inelastic), "%)")
    println(io, "│     ├── ionization        : ", b.ionization, "  (", pct(b.ionization), "%)")
    println(io, "│     ├── excitation        : ", b.excitation, "  (", pct(b.excitation), "%)")
    for (name, val) in b.inelastic_by_species
        println(io, "│     ├── ", rpad(name, 4), "             : ", val, "  (", pct(val), "%)")
    end
    println(io, "├── thermal heating         : ", b.heating, "  (", pct(b.heating), "%)")
    println(io, "├── escape (↑ top)          : ", b.escape, "  (", pct(b.escape),
            "%)  albedo = ", round(b.albedo, digits=3))
    println(io, "├── residual (subfloor+num) : ", b.residual, "  (", b.residual_fraction * 100, "%)")
    println(io, "└── ion-pairs = ", b.ionpairs, "   exc-events = ", b.excevents,
            "   z_centroid = ", round(b.z_centroid, digits=1), " km")
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
    z     = model.altitude_grid.h          # m (vertical altitude)
    s     = model.s_field                  # m (path length along the magnetic field line)
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
    # Vertical energy flux = Σ_beam Σ_E Ie·E·|μ| (matches field_aligned_beam_norm). The raw
    # (along-field, un-|μ|-weighted) Σ Ie·E is kept alongside for diagnostics.
    input  = 0.0; escape  = 0.0    # |μ|-weighted vertical flux (↓ / ↑)
    input_raw = 0.0; escape_raw = 0.0   # along-field, un-weighted
    @inbounds for iμ in 1:n_μ, iE in 1:n_E
        fe = Ie[n_z, iμ, iE] * E[iE]
        if μ[iμ] < 0
            input     += fe * abs(μ[iμ])
            input_raw += fe
        else
            escape     += fe * abs(μ[iμ])
            escape_raw += fe
        end
    end

    # ---- Omnidirectional flux for the volumetric (deposition / heating) terms -----------
    Ie_omni = dropdims(sum(Ie, dims=2), dims=2)            # [n_z, n_E]

    # ---- Inelastic energy deposition (∫ along the field line) ----------------------------
    # Σ_s Δs Σ_sp n Σ_levels threshold·Σ_E Ie_omni·σ, split into ionizing (≥1 secondary) and
    # non-ionizing channels, with the corresponding column reaction rates.
    dep_profile_total = zeros(n_z)      # energy-deposition profile [eV m⁻³ s⁻¹]
    ion_profile       = zeros(n_z)      # → ionization
    exc_profile       = zeros(n_z)      # → excitation
    ionpair_profile   = zeros(n_z)      # ion-pair production rate  [m⁻³ s⁻¹]
    excevent_profile  = zeros(n_z)      # excitation-event rate     [m⁻³ s⁻¹]
    inelastic_by_species = Pair{String,Float64}[]
    for sp in model.species
        σ      = sp.cross_sections          # [n_levels, n_E]
        levels = sp.excitation_levels       # [n_levels, 2]: (energy loss, #secondaries)
        dens   = sp.density                 # [n_z]
        dep_sp = zeros(n_z)
        for lvl in axes(levels, 1)
            E_loss = levels[lvl, 1]
            E_loss > 0 || continue          # skip elastic (no energy loss)
            secondaries = levels[lvl, 2]
            ionizing    = secondaries >= 1
            for iz in 1:n_z
                acc = 0.0
                @inbounds for iE in 1:n_E
                    acc += Ie_omni[iz, iE] * σ[lvl, iE]
                end
                rate = dens[iz] * acc       # reaction rate of this level [m⁻³ s⁻¹]
                dep  = rate * E_loss        # energy into this channel    [eV m⁻³ s⁻¹]
                dep_sp[iz] += dep
                if ionizing
                    ion_profile[iz]     += dep
                    ionpair_profile[iz] += rate * secondaries
                else
                    exc_profile[iz]      += dep
                    excevent_profile[iz] += rate
                end
            end
        end
        col = trapz_path(s, dep_sp)
        push!(inelastic_by_species, String(sp.name) => col)
        dep_profile_total .+= dep_sp
    end
    inelastic  = trapz_path(s, dep_profile_total)
    ionization = trapz_path(s, ion_profile)
    excitation = trapz_path(s, exc_profile)
    ionpairs   = trapz_path(s, ionpair_profile)
    excevents  = trapz_path(s, excevent_profile)

    # Energy-deposition centroid (vertical altitude, km): ∫ z·dep dz / ∫ dep dz.
    dep_col_z  = trapz_path(z, dep_profile_total)
    z_centroid = dep_col_z > 0 ? trapz_path(z, z .* dep_profile_total) / dep_col_z / 1e3 : NaN

    # ---- Thermal-electron heating (reuse the existing Coulomb-loss routine) --------------
    heating_profile = calculate_heating_rate(z, [0.0], reshape(Ie_omni, n_z, 1, n_E),
                                             E, ne, Te)[:, 1]
    heating = trapz_path(s, heating_profile)

    # ---- Closure ------------------------------------------------------------------------
    net      = input - escape
    residual = input - inelastic - heating - escape
    budget = EnergyBudget(input, escape, net, inelastic, ionization, excitation, heating,
                          residual, input > 0 ? residual / input : NaN,
                          input > 0 ? escape / input : NaN, input_raw, escape_raw,
                          ionpairs, excevents, z_centroid, inelastic_by_species)
    verbose && show(stdout, MIME"text/plain"(), budget)
    verbose && println()
    return budget
end

# Trapezoidal integral over a (possibly descending) 1-D grid x. Used both for the field-line
# path integral (x = model.s_field) and for the altitude integral of the deposition centroid.
function trapz_path(x, f)
    s = 0.0
    @inbounds for i in 1:(length(x) - 1)
        s += 0.5 * (f[i] + f[i + 1]) * abs(x[i + 1] - x[i])
    end
    return s
end
