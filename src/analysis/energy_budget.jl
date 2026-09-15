# ======================================================================================== #
#                              ENERGY-BUDGET DIAGNOSTIC                                     #
# ======================================================================================== #
#
# Global energy balance of the suprathermal electron population. At steady state the downward
# energy flux entering the top is accounted for by
#
#     E_in  =  E_inelastic  +  E_heating  +  E_escape  +  residual
#
# where
#   E_in        vertical energy flux carried by the *downward* beams at the top boundary
#               (= the imposed precipitation), Σ_beam Σ_E Ie·E·|μ|   [eV m⁻² s⁻¹]
#   E_escape    vertical energy flux carried by the *upward* beams at the top boundary
#               (backscattered electrons leaving the domain)
#   E_inelastic energy deposited in neutral excitation + ionization thresholds,
#               Σ_s w_s Σ_species n(s) Σ_levels threshold·Σ_E Ie_omni·σ
#   E_heating   energy transferred to the thermal electrons (Coulomb), via
#               `calculate_heating_rate`, integrated along the field line
#   residual    everything unaccounted: the energy of electrons degraded below the grid floor
#               (sub-floor thermalisation, a small *positive* term on a good grid) plus any
#               numerical non-conservation. A grid that violates the ΔE < threshold stability
#               bound shows the inelastic term ballooning and the residual going
#               large/negative (energy "created") — the failure this diagnostic catches.
#
# Volume integrals run ALONG THE MAGNETIC FIELD LINE (over `model.s_field = h / cos B_angle`),
# NOT over vertical altitude: the transport conserves the field-aligned flux (μ ∂Ie/∂s), so
# the column path length is s. Integrating over z instead injects a spurious 1/cos(B_angle)
# factor (it vanishes at B_angle = 0, where s == z).
#
# The quadrature is the solver's own (see `column_weights`), not a generic trapezoid rule, so
# the budget measures the energy the discrete operator actually deposits.
#
# The vertical-energy-flux projection (Σ Ie·E·|μ|) matches `field_aligned_beam_norm`, so
# E_in computed here equals the InputFlux's `IeE_tot` normalisation.
#
# A single time slice balances only at steady state. Over a time interval that starts and ends
# at rest the same identity holds for the time-integrated terms, which is what `trange` gives.

using Printf

"""
    EnergyBudget

Result of [`energy_budget`](@ref): where the precipitating energy flux ends up, with every
term column-integrated along the field line.

`interval` says which of the two forms this is. For a single time slice it is `nothing`, the
energy terms are fluxes in eV m⁻² s⁻¹ and `ionpairs` is a production rate in m⁻² s⁻¹. For a
time-integrated budget it is the `(t0, t1)` interval in seconds, the energy terms are
energies in eV m⁻² and `ionpairs` is a total count in m⁻².

Energy terms
- `input`, `escape`   |μ|-weighted *vertical* energy flux of the downward / upward beams at
  the top boundary (`input` matches `IeE_tot`).
- `net`               `input - escape`: the energy actually deposited in the column.
- `inelastic`         energy into neutral excitation + ionization thresholds.
- `ionization`, `excitation`  the inelastic term split into ionizing (≥1 secondary) and
  non-ionizing channels (`inelastic == ionization + excitation`).
- `heating`           energy transferred to the thermal electrons (Coulomb).
- `residual`          `input - inelastic - heating - escape` (sub-floor thermalisation +
  numerical non-conservation); `residual_fraction == residual / input`.
- `albedo`            `escape / input`.
- `inelastic_by_species`  the inelastic term per species, in model order.

Rates
- `ionpairs`          column ion-pair production.
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
    ionpairs::Float64
    inelastic_by_species::Vector{Pair{String,Float64}}
    interval::Union{Nothing,Tuple{Float64,Float64}}
end

# The scalar fields of EnergyBudget, in constructor order: the TOML writer and reader both
# iterate this tuple, so the file round-trips through the positional constructor.
const ENERGY_BUDGET_SCALAR_FIELDS = (
    :input, :escape, :net, :inelastic, :ionization, :excitation, :heating,
    :residual, :residual_fraction, :albedo, :ionpairs,
)

energy_units(b::EnergyBudget) = b.interval === nothing ? "eV m⁻² s⁻¹" : "eV m⁻²"
rate_units(b::EnergyBudget) = b.interval === nothing ? "m⁻² s⁻¹" : "m⁻²"

# Percentages, with two significant digits for a term below 0.1% so that a small but nonzero
# channel (thermal heating at a few keV) does not read as "0.0".
function percent_string(value, total)
    pct = 100 * value / total
    pct != 0 && abs(pct) < 0.1 && return string(round(pct; sigdigits = 2))
    return string(round(pct; digits = 1))
end

function Base.show(io::IO, ::MIME"text/plain", b::EnergyBudget)
    interval = b.interval
    span = interval === nothing ? "steady state" :
           "∫ over t = $(interval[1]) – $(interval[2]) s"
    value(x) = @sprintf("%.3g", x)
    row(name, x) = println(io, "  ", rpad(name, 21), lpad(value(x), 9), "  ",
                           lpad(percent_string(x, b.input), 9))

    println(io, "EnergyBudget — ", span, ", ∫ along the field line")
    println(io, rpad("input (↓ top)", 23), value(b.input), " ", energy_units(b))
    println(io, rpad("", 23), "value        % of input")
    row("ionization", b.ionization)
    row("excitation", b.excitation)
    row("thermal heating", b.heating)
    row("escape (↑ top)", b.escape)
    row("residual", b.residual)
    if !isempty(b.inelastic_by_species)
        shares = join(("$name $(percent_string(val, b.inelastic))"
                       for (name, val) in b.inelastic_by_species), ", ")
        println(io, "inelastic by species (% of inelastic): ", shares)
    end
    accounted = b.inelastic + b.heating + b.escape
    println(io, "albedo ", round(b.albedo; digits = 3),
            " · net energy per ion pair ", value(b.net / b.ionpairs), " eV",
            " · (deposited + backscattered)/input ", @sprintf("%.3f", accounted / b.input))
end

"""
    energy_budget(sim::AuroraSimulation; tidx=<last>, trange=nothing, verbose=true)
    energy_budget(sim_dir::AbstractString; tidx=nothing, trange=nothing, max_bytes=nothing,
                  verbose=true) -> EnergyBudget

Compute the energy balance of a finished simulation. Returns an [`EnergyBudget`](@ref); also
prints a summary unless `verbose=false`.

Two sources are accepted:
- a saved run directory `sim_dir`, reconstructing the model from
  `<sim_dir>/inputs/physics_state.jld2` (see [`load_model`](@ref)) and the electron flux from
  `<sim_dir>/simulation_data.nc`;
- an in-memory `sim`. For a steady-state run the flux comes straight from `sim.workspace.Ie`;
  for a time-dependent run the workspace holds only the last solver loop, so the call is
  forwarded to `sim.output.savedir`, which `run!` has written in full.

Pass at most one of `tidx` and `trange`:
- neither: the **last** time slice, which is the steady-state solution of an `n_t == 1` run;
- `tidx`: that single time slice;
- `trange`: a contiguous range of slices, or `:` for all of them, time-integrated with
  trapezoid weights over the time axis. The result then carries `interval = (t0, t1)` and
  holds energies (eV m⁻²) rather than fluxes.

A single slice balances only at steady state: on a time-dependent run energy is in transit,
and a *transient* input reads `input == 0` at any slice after the pulse has passed. A warning
points this out; `trange` is the form that closes for a transient that starts and ends at
rest, where `residual` also collects the energy left stored in the population.

Use as a guardrail: on a converged, stable run the residual is a small positive fraction (the
sub-floor thermalisation that AURORA does not track on-grid). A large or negative residual
flags energy non-conservation — e.g. a grid whose maximum bin width exceeds the lowest
ionization threshold, which destabilises the high→low energy-degradation sweep.

`max_bytes` caps the flux read from disk: one full time slice for a snapshot (no cap by
default, since the budget needs the whole slice) or one streaming chunk when integrating
(512 MiB by default; cf. [`foreach_Ie_time_chunk`](@ref)).
"""
function energy_budget(sim::AuroraSimulation; tidx = nothing, trange = nothing,
                       max_bytes::Union{Nothing,Real} = nothing, verbose::Bool = true)
    if sim.mode isa TimeDependentMode || trange !== nothing
        return energy_budget(sim.output.savedir; tidx, trange, max_bytes, verbose)
    end
    model  = sim.model
    Ie_raw = sim.workspace.Ie
    n_t    = size(Ie_raw, 2)
    it     = something(tidx, n_t)
    1 <= it <= n_t || throw(ArgumentError("tidx = $it out of range 1:$n_t"))
    warn_if_snapshot(n_t, it, verbose)
    n_z = length(model.altitude_grid.h)
    n_μ = length(model.pitch_angle_grid.μ_center)
    n_E = length(model.energy_grid.E_centers)
    # Reshape [n_z·n_μ, n_t, n_E] → [n_z, n_μ, n_E]; row = (i_μ-1)·n_z + i_z.
    Ie = reshape(@view(Ie_raw[:, it, :]), n_z, n_μ, n_E)
    return energy_budget_snapshot(model, Ie; verbose)
end

function energy_budget(sim_dir::AbstractString; tidx = nothing, trange = nothing,
                       max_bytes::Union{Nothing,Real} = nothing, verbose::Bool = true)
    tidx === nothing || trange === nothing ||
        throw(ArgumentError("pass either tidx (one time slice) or trange (a time " *
                            "integral), not both"))
    model = load_model(sim_dir)
    co = load_coordinates(sim_dir)
    trange === nothing &&
        return energy_budget_at(model, sim_dir, co, something(tidx, co.n_t),
                                something(max_bytes, Inf), verbose)
    return energy_budget_over(model, sim_dir, co, trange,
                              something(max_bytes, 512 * 1024^2), verbose)
end

# One time slice. A 300 keV slice is about 2.5 GiB, above `load_results`' generic 2 GiB safety
# default, and the budget needs the whole slice, so the cap is off unless the caller sets it.
function energy_budget_at(model, sim_dir, co, it::Integer, max_bytes, verbose)
    1 <= it <= co.n_t || throw(ArgumentError("tidx = $it out of range 1:$(co.n_t)"))
    warn_if_snapshot(co.n_t, it, verbose)
    res = load_results(sim_dir; tidx = it:it, max_bytes)
    Ie  = @view res.Ie[:, :, 1, :]                              # [n_z, n_μ, n_E]
    return energy_budget_snapshot(model, Ie; verbose)
end

# Trapezoidal time integral of the per-slice budgets over `trange`. The flux is streamed in
# time-chunks, so peak memory is bounded by `max_bytes` whatever the run length.
function energy_budget_over(model, sim_dir, co, trange, max_bytes, verbose)
    ts = trange === Colon() ? (1:co.n_t) : trange
    (ts isa AbstractUnitRange{<:Integer} && first(ts) >= 1 && last(ts) <= co.n_t) ||
        throw(ArgumentError("trange must be a Colon (:) or a contiguous range within " *
                            "1:$(co.n_t)"))
    length(ts) >= 2 ||
        throw(ArgumentError("need ≥ 2 time slices to integrate; pass tidx for a single " *
                            "snapshot"))
    t = co.t[ts]
    w = trapz_weights(t)               # ∫ f dt ≈ Σ w[k] f[k]

    names = [String(sp.name) for sp in model.species]
    sums  = Dict(n => 0.0 for n in names)
    input = 0.0; escape = 0.0; inelastic = 0.0; ionization = 0.0
    excitation = 0.0; heating = 0.0; ionpairs = 0.0
    # `foreach_Ie_time_chunk` reuses its buffer between calls, so each slice is reduced to a
    # budget before the next chunk is read.
    foreach_Ie_time_chunk(sim_dir; trange = ts, max_bytes) do Ie_chunk, t_range
        for (j, it) in enumerate(t_range)
            b = energy_budget_snapshot(model, @view(Ie_chunk[:, :, j, :]); verbose = false)
            wk = w[it - first(ts) + 1]
            input += wk * b.input;           escape += wk * b.escape
            inelastic += wk * b.inelastic;   ionization += wk * b.ionization
            excitation += wk * b.excitation; heating += wk * b.heating
            ionpairs += wk * b.ionpairs
            for (name, val) in b.inelastic_by_species
                sums[name] += wk * val
            end
        end
    end

    budget = close_budget(input, escape, inelastic, ionization, excitation, heating,
                          ionpairs, [n => sums[n] for n in names],
                          (first(t), last(t)))
    verbose && (show(stdout, MIME"text/plain"(), budget); println())
    return budget
end

"""
    make_energy_budget_file(sim_or_dir; tidx=nothing, trange=nothing, verbose=true)
        -> EnergyBudget

Compute the energy budget with [`energy_budget`](@ref) (same keywords) and save it to
`<savedir>/analysis/energy_budget.toml`. The compact TOML file holds every scalar field of
[`EnergyBudget`](@ref), the species-resolved inelastic terms in model order, the time
interval, and the units, so the budget stays readable when the much larger
`simulation_data.nc` is moved or deleted. Read it back with [`load_energy_budget`](@ref).
"""
function make_energy_budget_file(sim::AuroraSimulation; kwargs...)
    budget = energy_budget(sim; kwargs...)
    return write_energy_budget_file(sim.output.savedir, budget)
end

function make_energy_budget_file(sim_dir::AbstractString; kwargs...)
    budget = energy_budget(sim_dir; kwargs...)
    return write_energy_budget_file(sim_dir, budget)
end

function write_energy_budget_file(sim_dir::AbstractString, budget::EnergyBudget)
    analysis_dir = joinpath(sim_dir, "analysis")
    mkpath(analysis_dir)
    savefile = joinpath(analysis_dir, "energy_budget.toml")
    values = Dict(String(name) => getfield(budget, name) for name in ENERGY_BUDGET_SCALAR_FIELDS)
    data = Dict{String,Any}(
        "schema_version" => 1,
        "energy_units" => energy_units(budget),
        "rate_units" => rate_units(budget),
        "values" => values,
        # An array of tables, so the species keep their model order on the way back in.
        "inelastic_by_species" => [Dict("species" => name, "value" => val)
                                   for (name, val) in budget.inelastic_by_species],
    )
    budget.interval === nothing || (data["interval"] = collect(budget.interval))

    # Write through a temporary file in the same directory, so an interrupted write leaves
    # any previous energy_budget.toml intact rather than truncated.
    tmp, io = mktemp(analysis_dir)
    try
        TOML.print(io, data; sorted = true)
        close(io)
        mv(tmp, savefile; force = true)
    finally
        isopen(io) && close(io)
        isfile(tmp) && rm(tmp; force = true)
    end
    println("Energy budget saved in $savefile")
    return budget
end

"""
    load_energy_budget(sim_or_directory) -> EnergyBudget

Load the compact result written by [`make_energy_budget_file`](@ref). Accepts either a saved
run directory or an `AuroraSimulation` whose output directory contains the file.
"""
function load_energy_budget(sim_dir::AbstractString)
    path = joinpath(sim_dir, "analysis", "energy_budget.toml")
    isfile(path) || throw(ArgumentError("no energy_budget.toml found under $(dirname(path))"))
    check_analysis_freshness(path, sim_dir)
    data = TOML.parsefile(path)
    get(data, "schema_version", nothing) == 1 ||
        throw(ArgumentError("unsupported energy-budget schema in $path"))
    values = data["values"]
    scalars = (Float64(values[String(name)]) for name in ENERGY_BUDGET_SCALAR_FIELDS)
    entries = data["inelastic_by_species"]
    entries isa AbstractVector ||
        throw(ArgumentError("$path stores inelastic_by_species as $(typeof(entries)); it " *
                            "must be an array of {species, value} tables. Rewrite the file " *
                            "with make_energy_budget_file."))
    species = [String(entry["species"]) => Float64(entry["value"]) for entry in entries]
    interval = haskey(data, "interval") ?
               (Float64(data["interval"][1]), Float64(data["interval"][2])) : nothing
    return EnergyBudget(scalars..., species, interval)
end

load_energy_budget(sim::AuroraSimulation) = load_energy_budget(sim.output.savedir)

# A single slice balances only at steady state. Warn (when printing) that a snapshot of a
# time-dependent run does not conserve, since a transient input can read input == 0 at a slice
# taken after the pulse — the usual source of a "surprising" zero/blown-up budget.
function warn_if_snapshot(n_t, tidx, verbose)
    if n_t > 1 && verbose
        @warn "energy_budget is a single-time snapshot (tidx = $tidx of $n_t); for a " *
              "time-dependent run the balance does not close (energy in transit), and a " *
              "transient input reads 0 after the pulse. Pass `tidx` to pick a time, or " *
              "`trange` to integrate over time."
    end
end

# Assemble the derived terms shared by the snapshot and the time integral.
function close_budget(input, escape, inelastic, ionization, excitation, heating, ionpairs,
                      inelastic_by_species, interval)
    net      = input - escape
    residual = input - inelastic - heating - escape
    return EnergyBudget(input, escape, net, inelastic, ionization, excitation, heating,
                        residual, input > 0 ? residual / input : NaN,
                        input > 0 ? escape / input : NaN, ionpairs,
                        inelastic_by_species, interval)
end

# Budget of one flux snapshot `Ie`, shape [n_z, n_μ, n_E].
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
    # The loops below index z, μ, E and Ie with the same indices, so the axes must agree.
    axes(Ie) == (eachindex(z), eachindex(μ), eachindex(E)) ||
        throw(ArgumentError("Ie snapshot $(size(Ie)) does not match the model grid " *
                            "(n_z, n_μ, n_E) = ($n_z, $n_μ, $n_E)"))

    # ---- Boundary energy fluxes at the top altitude --------------------------------------
    # Vertical energy flux = Σ_beam Σ_E Ie·E·|μ| (matches field_aligned_beam_norm).
    input = 0.0; escape = 0.0
    i_top = lastindex(Ie, 1)
    for iμ in axes(Ie, 2), iE in axes(Ie, 3)
        fe = Ie[i_top, iμ, iE] * E[iE] * abs(μ[iμ])
        if μ[iμ] < 0
            input += fe
        else
            escape += fe
        end
    end

    # ---- Weighted omnidirectional flux for the volumetric (deposition / heating) terms ---
    # Every volumetric term is linear in the local flux, so folding each beam's quadrature
    # weight into the flux integrates all of them in one pass: the profiles below are already
    # weighted, and their plain sums over altitude are the column integrals.
    w = column_weights(s)
    Ie_omni = zeros(n_z, n_E)                              # [n_z, n_E], weight-folded
    for iμ in axes(Ie, 2)
        w_μ = μ[iμ] < 0 ? w.down : w.up
        for iE in axes(Ie, 3), iz in axes(Ie, 1)
            Ie_omni[iz, iE] += w_μ[iz] * Ie[iz, iμ, iE]
        end
    end

    # ---- Inelastic energy deposition (∫ along the field line) ----------------------------
    # Σ_s w_s Σ_sp n Σ_levels threshold·Σ_E Ie_omni·σ, split into ionizing (≥1 secondary) and
    # non-ionizing channels, with the ion-pair production rate alongside.
    dep_profile_total = zeros(n_z)      # weighted energy deposition [eV m⁻² s⁻¹]
    ion_profile       = zeros(n_z)      # → ionization
    exc_profile       = zeros(n_z)      # → excitation
    ionpair_profile   = zeros(n_z)      # weighted ion-pair production [m⁻² s⁻¹]
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
            for iz in eachindex(dep_sp, dens)
                acc = 0.0
                for iE in eachindex(E)
                    acc += Ie_omni[iz, iE] * σ[lvl, iE]
                end
                rate = dens[iz] * acc       # reaction rate of this level [m⁻³ s⁻¹]
                dep  = rate * E_loss        # energy into this channel    [eV m⁻³ s⁻¹]
                dep_sp[iz] += dep
                if ionizing
                    ion_profile[iz]     += dep
                    ionpair_profile[iz] += rate * secondaries
                else
                    exc_profile[iz] += dep
                end
            end
        end
        push!(inelastic_by_species, String(sp.name) => sum(dep_sp))
        dep_profile_total .+= dep_sp
    end
    inelastic  = sum(dep_profile_total)
    ionization = sum(ion_profile)
    excitation = sum(exc_profile)
    ionpairs   = sum(ionpair_profile)

    # ---- Thermal-electron heating (reuse the existing Coulomb-loss routine) --------------
    heating_profile = calculate_heating_rate(z, [0.0], reshape(Ie_omni, n_z, 1, n_E),
                                             E, ne, Te)[:, 1]
    heating = sum(heating_profile)

    budget = close_budget(input, escape, inelastic, ionization, excitation, heating,
                          ionpairs, inelastic_by_species, nothing)
    verbose && (show(stdout, MIME"text/plain"(), budget); println())
    return budget
end

# Quadrature weights for the column integral along the field line, taken from the solvers so
# that summing the discrete equations with them telescopes the transport term into the
# boundary fluxes and leaves exactly the sinks this budget measures.
#
# Both solvers difference μ ∂Ie/∂s upwind, with the interval the beam comes from in the
# denominator (`build_spatial_operators` in src/solvers/sparse_indexing.jl): downward beams
# use s[i+1] - s[i], upward beams s[i] - s[i-1]. The collision sinks sit in the same rows
# with no cell length of their own, so they carry the same weight. The first and last rows
# hold boundary conditions rather than a balance, so they get no weight: the flux at the top
# is already counted as `input`/`escape`, and the flux the bottom row zeroes leaves the
# domain without a budget term.
function column_weights(s)
    down = zeros(length(s))
    up   = zeros(length(s))
    lo, hi = firstindex(s), lastindex(s)
    for k in (lo + 1):(hi - 1)
        down[k - lo + 1] = abs(s[k + 1] - s[k])
        up[k - lo + 1]   = abs(s[k] - s[k - 1])
    end
    return (; down, up)
end

# Trapezoidal quadrature weights w for samples at points x, so that ∫ f dx ≈ Σ w[k] f[k].
function trapz_weights(x)
    w = zeros(length(x))
    length(x) == 1 && return w
    lo_idx, hi_idx = firstindex(x), lastindex(x)
    for (j, k) in enumerate(eachindex(x))
        lo = k == lo_idx ? x[lo_idx] : x[k - 1]
        hi = k == hi_idx ? x[hi_idx] : x[k + 1]
        w[j] = (hi - lo) / 2
    end
    return w
end
