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
#
# Works either on an in-memory `sim` or on a saved run directory (reconstructing the model and
# flux from disk) — see the two `energy_budget` methods below.

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
    TimeIntegratedEnergyBudget

Result of [`energy_budget_integrated`](@ref): an [`EnergyBudget`](@ref) whose terms have been
integrated over a time interval `[t0, t1]` (seconds). Energy fluxes are therefore in **eV m⁻²**
(energy per area over the interval) and `ionpairs`/`excevents` are **total counts m⁻²**. All
`EnergyBudget` fields are forwarded (e.g. `b.input`, `b.residual`); additionally `t0`, `t1`, and
`span` (= `t1 - t0`) are available.

For a transient that starts and ends at rest the balance closes, and `residual` then collects
the energy left stored in the population over the interval plus numerical non-conservation.
"""
struct TimeIntegratedEnergyBudget
    budget::EnergyBudget
    t0::Float64
    t1::Float64
end

function Base.getproperty(b::TimeIntegratedEnergyBudget, s::Symbol)
    s === :span && return getfield(b, :t1) - getfield(b, :t0)
    (s === :budget || s === :t0 || s === :t1) && return getfield(b, s)
    return getproperty(getfield(b, :budget), s)   # forward EnergyBudget fields
end

Base.propertynames(b::TimeIntegratedEnergyBudget) =
    (:t0, :t1, :span, propertynames(getfield(b, :budget))...)

function Base.show(io::IO, ::MIME"text/plain", b::TimeIntegratedEnergyBudget)
    pct(x) = round(100x / b.input, digits=2)
    println(io, "TimeIntegratedEnergyBudget (eV m⁻², ∫ over Δt = ", round(b.span, digits=4),
            " s, ∫ along field line):")
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
    println(io, "└── residual (Δstored+num)  : ", b.residual, "  (", b.residual_fraction * 100, "%)")
end

"""
    energy_budget(sim::AuroraSimulation; tidx=<last>, verbose=true) -> EnergyBudget
    energy_budget(sim_dir::AbstractString; tidx=<last>, verbose=true) -> EnergyBudget

Compute the steady-state energy balance of a finished simulation. Returns an
[`EnergyBudget`](@ref); also prints a summary unless `verbose=false`.

Two sources are accepted:
- an in-memory `sim`, read from `sim.cache.Ie`;
- a saved run directory `sim_dir`, reconstructing the model from
  `<sim_dir>/inputs/physics_state.jld2` (see [`load_model`](@ref)) and the electron flux from
  `<sim_dir>/simulation_data.nc`.

The budget is evaluated at a single time slice, `tidx` (default: the **last** slice, which is
the steady-state solution for an `n_t == 1` run). This balance closes only at steady state: on
a **time-dependent** run a single snapshot does not conserve (energy is in transit / stored),
and in particular a *transient* input (e.g. a pulse) will read `input == 0` at any time after
the pulse has passed — pass `tidx` to evaluate at the time of interest (e.g. the input peak). A
warning is emitted when the run has more than one time step.

Use as a guardrail: on a converged, stable run the residual is a small positive fraction (the
sub-floor thermalisation that AURORA does not track on-grid). A large or negative residual
flags energy non-conservation — e.g. a grid whose maximum bin width exceeds the lowest
ionization threshold, which destabilises the high→low energy-degradation sweep.
"""
function energy_budget(sim::AuroraSimulation; tidx::Integer = size(sim.cache.Ie, 2),
                       verbose::Bool = true)
    model  = sim.model
    Ie_raw = sim.cache.Ie
    n_t    = size(Ie_raw, 2)
    1 <= tidx <= n_t || throw(ArgumentError("tidx = $tidx out of range 1:$n_t"))
    warn_if_time_dependent(n_t, tidx, verbose)
    n_z = length(model.altitude_grid.h)
    n_μ = length(model.pitch_angle_grid.μ_center)
    n_E = length(model.energy_grid.E_centers)
    # Reshape [n_z·n_μ, n_t, n_E] → [n_z, n_μ, n_E]; row = (i_μ-1)·n_z + i_z.
    Ie = reshape(@view(Ie_raw[:, tidx, :]), n_z, n_μ, n_E)
    return energy_budget_snapshot(model, Ie; verbose)
end

function energy_budget(sim_dir::AbstractString; tidx::Union{Nothing,Integer} = nothing,
                       verbose::Bool = true)
    model = load_model(sim_dir)
    # Recompute the field-line path length from the saved geometry, so the diagnostic is robust
    # to older physics_state.jld2 files written before `s_field` existed.
    model.s_field = model.altitude_grid.h ./ cosd(model.B_angle_to_zenith)
    n_t = load_coordinates(sim_dir).n_t
    it  = something(tidx, n_t)              # default: final slice
    1 <= it <= n_t || throw(ArgumentError("tidx = $it out of range 1:$n_t"))
    warn_if_time_dependent(n_t, it, verbose)
    res = load_results(sim_dir; tidx = it:it)
    Ie  = @view res.Ie[:, :, 1, :]                              # [n_z, n_μ, n_E]
    return energy_budget_snapshot(model, Ie; verbose)
end

"""
    energy_budget_integrated(sim_dir; trange=:, max_bytes=512*1024^2, verbose=true)
        -> TimeIntegratedEnergyBudget

Time-integrated energy balance for a *time-dependent* run, computed from a saved directory by
trapezoidally integrating the per-time-slice [`energy_budget`](@ref) over the time axis (or the
contiguous integer sub-range `trange`). Returns a [`TimeIntegratedEnergyBudget`](@ref) (energy
fluxes in eV m⁻²; rates as total counts m⁻²); also prints a summary unless `verbose=false`.

Unlike the single-slice `energy_budget`, this closes for a transient that starts and ends at
rest: `input = inelastic + heating + escape + residual`, with `residual` then collecting the
energy left stored in the population over `[t0, t1]` plus numerical non-conservation. This is the
appropriate tool for a pulsed / Alfvénic precipitation run, where a single snapshot does not
conserve.

The flux is streamed from `simulation_data.nc` in time-chunks (a single pass, peak memory
bounded by `max_bytes`; cf. [`foreach_Ie_time_chunk`](@ref)) rather than loaded all at once. A
sub-range `trange` still streams the whole file and filters in memory. `z_centroid` is the
deposition-energy-weighted average of the per-slice centroids.
"""
function energy_budget_integrated(sim_dir::AbstractString; trange = Colon(),
                                  max_bytes::Real = 512 * 1024^2, verbose::Bool = true)
    model = load_model(sim_dir)
    model.s_field = model.altitude_grid.h ./ cosd(model.B_angle_to_zenith)
    co = load_coordinates(sim_dir)
    ts = trange === Colon() ? (1:co.n_t) : trange
    (ts isa AbstractUnitRange{<:Integer} && first(ts) >= 1 && last(ts) <= co.n_t) ||
        throw(ArgumentError("trange must be a Colon (:) or a contiguous range within 1:$(co.n_t)"))
    length(ts) >= 2 ||
        throw(ArgumentError("need ≥ 2 time slices to integrate; use `energy_budget` for a single snapshot"))
    tt = co.t[ts]
    w  = trapz_weights(tt)             # ∫ f dt ≈ Σ w[k] f[k]
    lo, hi = first(ts), last(ts)

    # Stream the flux in time-chunks (a single pass over simulation_data.nc, peak memory bounded
    # by `max_bytes`) and keep the cheap per-slice snapshot budgets; integrate them afterwards.
    # `push!` mutates the vectors in place, so the streaming closure boxes none of the totals.
    budgets = EnergyBudget[]
    idxs    = Int[]
    foreach_Ie_time_chunk(sim_dir; max_bytes) do Ie_chunk, t_range
        for (j, it) in enumerate(t_range)
            lo <= it <= hi || continue
            push!(budgets, energy_budget_snapshot(model, @view(Ie_chunk[:, :, j, :]);
                                                  verbose = false))
            push!(idxs, it)
        end
    end

    names = [String(sp.name) for sp in model.species]
    sums  = Dict(n => 0.0 for n in names)
    input = 0.0; escape = 0.0; inelastic = 0.0; ionization = 0.0; excitation = 0.0
    heating = 0.0; input_raw = 0.0; escape_raw = 0.0; ionpairs = 0.0; excevents = 0.0
    cz_num = 0.0; cz_den = 0.0
    for (m, b) in enumerate(budgets)
        wk = w[idxs[m] - lo + 1]
        input += wk*b.input;          escape += wk*b.escape
        inelastic += wk*b.inelastic;  ionization += wk*b.ionization
        excitation += wk*b.excitation; heating += wk*b.heating
        input_raw += wk*b.input_raw;  escape_raw += wk*b.escape_raw
        ionpairs += wk*b.ionpairs;    excevents += wk*b.excevents
        for (name, val) in b.inelastic_by_species
            sums[name] += wk*val
        end
        if b.inelastic > 0 && isfinite(b.z_centroid)
            cz_num += wk*b.inelastic*b.z_centroid; cz_den += wk*b.inelastic
        end
    end
    net        = input - escape
    residual   = input - inelastic - heating - escape
    z_centroid = cz_den > 0 ? cz_num / cz_den : NaN
    budget = EnergyBudget(input, escape, net, inelastic, ionization, excitation, heating,
                          residual, input > 0 ? residual / input : NaN,
                          input > 0 ? escape / input : NaN, input_raw, escape_raw,
                          ionpairs, excevents, z_centroid, [n => sums[n] for n in names])
    result = TimeIntegratedEnergyBudget(budget, first(tt), last(tt))
    verbose && (show(stdout, MIME"text/plain"(), result); println())
    return result
end

# The budget closes only at steady state. Warn (once per call, when printing) that a snapshot of
# a time-dependent run does not conserve, since a transient input can read input == 0 at a slice
# taken after the pulse — the usual source of a "surprising" zero/blown-up budget.
function warn_if_time_dependent(n_t, tidx, verbose)
    if n_t > 1 && verbose
        @warn "energy_budget is a single-time snapshot (tidx = $tidx of $n_t); for a " *
              "time-dependent run the balance does not close (energy in transit), and a " *
              "transient input reads 0 after the pulse. Pass `tidx` to pick a time."
    end
end

# Core computation shared by both `energy_budget` methods. `Ie` is the steady-state flux
# snapshot, shape [n_z, n_μ, n_E].
function energy_budget_snapshot(model, Ie; verbose::Bool = true)
    eg  = model.energy_grid
    z   = model.altitude_grid.h          # m (vertical altitude)
    s   = model.s_field                  # m (path length along the magnetic field line)
    μ   = model.pitch_angle_grid.μ_center
    E   = eg.E_centers                   # eV
    ne  = model.ionosphere.ne
    Te  = model.ionosphere.Te
    n_z = length(z)
    n_μ = length(μ)
    n_E = length(E)
    size(Ie) == (n_z, n_μ, n_E) ||
        throw(ArgumentError("Ie snapshot $(size(Ie)) does not match the model grid " *
                            "(n_z, n_μ, n_E) = ($n_z, $n_μ, $n_E)"))

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

# Trapezoidal quadrature weights w for samples at points x, so that ∫ f dx ≈ Σ w[k] f[k].
# Lets the time-integrated budget integrate each per-slice quantity with one weight vector.
function trapz_weights(x)
    n = length(x)
    w = zeros(n)
    n == 1 && return w
    @inbounds for k in 1:n
        lo = k == 1 ? x[1] : x[k - 1]
        hi = k == n ? x[n] : x[k + 1]
        w[k] = (hi - lo) / 2
    end
    return w
end
