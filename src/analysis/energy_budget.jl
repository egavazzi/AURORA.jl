# ======================================================================================== #
#                              ENERGY-BUDGET DIAGNOSTIC                                     #
# ======================================================================================== #
#
# Global energy balance of the suprathermal electron population. At steady state the downward
# energy flux entering the top is accounted for by
#
#     E_in  =  E_inelastic  +  E_heating  +  E_subfloor  +  E_escape  +  E_bottom  +  residual
#
# (the terms are defined in the `EnergyBudget` docstring). A single time slice balances only
# at steady state. Over a time interval that starts and ends at rest the same identity holds
# for the time-integrated terms, which is what `trange` gives.

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
- `bottom_escape`     the downward energy flux escaping at the bottom of the grid.
- `net`               `input - escape - bottom_escape`: the energy deposited in the modelled
  column, and the numerator of the energy per ion pair.
- `inelastic`         energy into neutral excitation + ionization thresholds.
- `ionization`, `excitation`  the inelastic term split into ionizing (≥1 secondary) and
  non-ionizing channels (`inelastic == ionization + excitation`).
- `heating`           energy transferred to the thermal electrons (Coulomb).
- `subfloor`          energy of the electrons that inelastic collisions and the Coulomb loss
  take below the lowest edge of the energy grid, where they thermalise.
- `residual`          `input - inelastic - heating - subfloor - escape - bottom_escape`,
  the numerical non-conservation of the scheme.
  `residual_fraction == residual / input`.
- `albedo`            `escape / input`.
- `inelastic_by_species`  the inelastic term per species, in model order.

Rates
- `ionpairs`          column ion-pair production.
"""
struct EnergyBudget
    input::Float64
    escape::Float64
    bottom_escape::Float64
    net::Float64
    inelastic::Float64
    ionization::Float64
    excitation::Float64
    heating::Float64
    subfloor::Float64
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
    :input, :escape, :bottom_escape, :net, :inelastic, :ionization, :excitation, :heating,
    :subfloor, :residual, :residual_fraction, :albedo, :ionpairs,
)

energy_units(b::EnergyBudget) = b.interval === nothing ? "eV m⁻² s⁻¹" : "eV m⁻²"
rate_units(b::EnergyBudget) = b.interval === nothing ? "m⁻² s⁻¹" : "m⁻²"

# Percentages, with two significant digits for a term below 0.1% so that a small but nonzero
# channel (thermal heating at a few keV) does not read as "0.0".
function percent_string(value, total)
    pct = 100 * value / total
    # @sprintf rather than round: round(x; sigdigits = 2) can sometimes return a float whose
    # shortest representation is not two digits, e.g. 8.51e-24 rounds to 8.499999999999999e-24.
    pct != 0 && abs(pct) < 0.1 && return replace(@sprintf("%.2g", pct), "e-0" => "e-")
    return string(round(pct; digits = 1))
end

function Base.show(io::IO, ::MIME"text/plain", b::EnergyBudget)
    interval = b.interval
    span = interval === nothing ? "steady state" :
           "∫ over t = $(interval[1]) – $(interval[2]) s"
    value(x) = replace(@sprintf("%.3g", x), "e+" => "e")   # 2.95e16, not 2.95e+16
    row(name, x) = println(io, "  ", rpad(name, 24), lpad(value(x), 9), "  ",
                           lpad(percent_string(x, b.input), 9))

    println(io, "EnergyBudget — ", span, ", ∫ along the field line")
    println(io, rpad("input (↓ top)", 26), value(b.input), " ", energy_units(b))
    println(io, rpad("", 26), "value        % of input")
    row("ionization", b.ionization)
    row("excitation", b.excitation)
    row("thermal heating", b.heating)
    row("sub-floor thermalisation", b.subfloor)
    row("escape (↑ top)", b.escape)
    row("escape (↓ bottom)", b.bottom_escape)
    row("residual", b.residual)
    if !isempty(b.inelastic_by_species)
        shares = join(("$name $(percent_string(val, b.inelastic))"
                       for (name, val) in b.inelastic_by_species), ", ")
        println(io, "inelastic by species (% of inelastic): ", shares)
    end
    accounted = b.inelastic + b.heating + b.subfloor + b.escape + b.bottom_escape
    println(io, "albedo ", round(b.albedo; digits = 3),
            " · net energy per ion pair ", value(b.net / b.ionpairs), " eV",
            " · (deposited + escaped)/input ", @sprintf("%.3f", accounted / b.input))
end

"""
    energy_budget(sim::AuroraSimulation; tidx=<last>, t=nothing, trange=nothing, verbose=true)
    energy_budget(sim_dir::AbstractString; tidx=nothing, t=nothing, trange=nothing,
                  max_bytes=nothing, verbose=true) -> EnergyBudget

Compute the energy balance of a finished simulation. Returns an [`EnergyBudget`](@ref); also
prints a summary unless `verbose=false`.

Accepts a saved run directory or an in-memory simulation.

Pass at most one of `tidx`, `t` and `trange`:
- none: the last saved slice, which is the steady state of a steady-state run;
- `tidx`: a single time slice;
- `t`: the single saved slice nearest that time in seconds (the earlier one on a tie);
- `trange`: several slices, time-integrated. It can take `:` for every saved slice, a
    contiguous range of slice indices, or a `(t0, t1)` tuple of times in seconds which
    selects every slice with `t0 <= t <= t1`. The result then holds energies (eV m⁻²)
    rather than fluxes.
"""
function energy_budget(sim::AuroraSimulation; tidx = nothing, t = nothing, trange = nothing,
                       max_bytes::Union{Nothing,Real} = nothing, verbose::Bool = true)
    if sim.mode isa TimeDependentMode || trange !== nothing || t !== nothing
        return energy_budget(sim.output.savedir; tidx, t, trange, max_bytes, verbose)
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

function energy_budget(sim_dir::AbstractString; tidx = nothing, t = nothing,
                       trange = nothing, max_bytes::Union{Nothing,Real} = nothing,
                       verbose::Bool = true)
    count(!isnothing, (tidx, t, trange)) <= 1 ||
        throw(ArgumentError("pass at most one of tidx (one slice by index), t (one slice " *
                            "by time) and trange (a time integral)"))
    model = load_model(sim_dir)
    co = load_coordinates(sim_dir)
    t === nothing || (tidx = time_index_nearest(co, t))
    if trange === nothing
        # A slice must be read whole, and a 300 keV slice (about 2.5 GiB) exceeds the 2 GiB
        # default cap of `load_results`, so the cap is off unless the caller set one.
        return energy_budget_at(model, sim_dir, co, something(tidx, co.n_t),
                                something(max_bytes, Inf), verbose)
    else
        return energy_budget_over(model, sim_dir, co, trange,
                              something(max_bytes, 512 * 1024^2), verbose)
    end
end


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
    ts = time_index_range(co, trange)
    length(ts) >= 2 ||
        throw(ArgumentError("need ≥ 2 time slices to integrate; pass tidx for a single " *
                            "snapshot"))
    t = co.t[ts]
    w = trapz_weights(t)               # ∫ f dt ≈ Σ w[k] f[k]

    names = [String(sp.name) for sp in model.species]
    sums  = Dict(n => 0.0 for n in names)
    input = 0.0; escape = 0.0; bottom_escape = 0.0; inelastic = 0.0; ionization = 0.0
    excitation = 0.0; heating = 0.0; subfloor = 0.0; ionpairs = 0.0
    subfloor_factors = budget_subfloor_factors(model)
    # `foreach_Ie_time_chunk` reuses its buffer between calls, so each slice is reduced to a
    # budget before the next chunk is read.
    foreach_Ie_time_chunk(sim_dir; trange = ts, max_bytes) do Ie_chunk, t_range
        for (j, it) in enumerate(t_range)
            b = energy_budget_snapshot(model, @view(Ie_chunk[:, :, j, :]); verbose = false,
                                       subfloor_factors)
            wk = w[it - first(ts) + 1]
            input += wk * b.input;           escape += wk * b.escape
            bottom_escape += wk * b.bottom_escape
            inelastic += wk * b.inelastic;   ionization += wk * b.ionization
            excitation += wk * b.excitation; heating += wk * b.heating
            subfloor += wk * b.subfloor;     ionpairs += wk * b.ionpairs
            for (name, val) in b.inelastic_by_species
                sums[name] += wk * val
            end
        end
    end

    budget = close_budget(input, escape, bottom_escape, inelastic, ionization, excitation,
                          heating, subfloor, ionpairs, [n => sums[n] for n in names],
                          (first(t), last(t)))
    verbose && (show(stdout, MIME"text/plain"(), budget); println())
    return budget
end

"""
    make_energy_budget_file(sim_or_dir; tidx=nothing, trange=nothing, verbose=true)
        -> EnergyBudget

Compute the energy budget with [`energy_budget`](@ref) (same keywords) and save it to
`<savedir>/analysis/energy_budget.toml`. The TOML file holds every scalar field of
[`EnergyBudget`](@ref), the species-resolved inelastic terms in model order, the time
interval, and the units. Read it back with [`load_energy_budget`](@ref).
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

Load the result written by [`make_energy_budget_file`](@ref). Accepts either a saved
run directory or an `AuroraSimulation` whose output directory contains the file.
"""
function load_energy_budget(sim_dir::AbstractString)
    path = joinpath(sim_dir, "analysis", "energy_budget.toml")
    isfile(path) || throw(ArgumentError("no energy_budget.toml found under $(dirname(path))"))
    check_analysis_freshness(path, sim_dir)
    data = TOML.parsefile(path)
    values = data["values"]
    missing_fields = [String(name) for name in ENERGY_BUDGET_SCALAR_FIELDS
                      if !haskey(values, String(name))]
    isempty(missing_fields) ||
        throw(ArgumentError("$path has no value for $(join(missing_fields, ", ")), so it " *
                            "does not describe the current EnergyBudget. Rewrite the file " *
                            "with make_energy_budget_file."))
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

# The saved slice closest in time to `t`, the earlier one when two are equally close.
function time_index_nearest(co, t)
    issorted(co.t) ||
        throw(ArgumentError("the saved time axis is not sorted, so there is no nearest " *
                            "slice to a given time"))
    first(co.t) <= t <= last(co.t) ||
        throw(ArgumentError("t = $t lies outside the saved time span " *
                            "$(first(co.t)) – $(last(co.t)) s"))
    return argmin(abs.(co.t .- t))
end

# Resolve the `trange` selector into the contiguous range of saved time indices it names:
# `:` is every slice, an integer range is itself, and a `(t0, t1)` tuple is every slice whose
# time lies in the closed interval, which is contiguous because the time axis is sorted.
function time_index_range(co, trange)
    trange === Colon() && return 1:co.n_t
    if trange isa AbstractUnitRange{<:Integer}
        first(trange) >= 1 && last(trange) <= co.n_t ||
            throw(ArgumentError("trange = $trange is out of bounds for the $(co.n_t) saved " *
                                "time slices"))
        return trange
    end
    if trange isa Tuple{Real,Real}
        t0, t1 = trange
        t0 <= t1 ||
            throw(ArgumentError("trange = $trange runs backwards; pass (t0, t1) with " *
                                "t0 <= t1"))
        issorted(co.t) ||
            throw(ArgumentError("the saved time axis is not sorted, so a time interval does " *
                                "not select a contiguous range of slices"))
        inside = findall(t -> t0 <= t <= t1, co.t)
        length(inside) >= 2 ||
            throw(ArgumentError("trange = $trange selects $(length(inside)) of the " *
                                "$(co.n_t) saved time slices, which span " *
                                "$(first(co.t)) – $(last(co.t)) s; at least 2 are needed " *
                                "to integrate"))
        return first(inside):last(inside)
    end
    throw(ArgumentError("trange must be a Colon (:), a contiguous range of time indices, " *
                        "or a (t0, t1) tuple of times in seconds; got $(typeof(trange))"))
end

# A single slice balances only at steady state. Warn (when printing) that a snapshot of a
# time-dependent run does not conserve.
function warn_if_snapshot(n_t, tidx, verbose)
    if n_t > 1 && verbose
        @warn "energy_budget is a single-time snapshot (tidx = $tidx of $n_t); for a " *
              "time-dependent run the balance does not close (energy in transit), and a " *
              "transient input reads 0 after the pulse. Pass `tidx` to pick a time, or " *
              "`trange` to integrate over time."
    end
end

# Assemble the derived terms shared by the snapshot and the time integral.
function close_budget(input, escape, bottom_escape, inelastic, ionization, excitation,
                      heating, subfloor, ionpairs, inelastic_by_species, interval)
    net      = input - escape - bottom_escape
    residual = input - inelastic - heating - subfloor - escape - bottom_escape
    return EnergyBudget(input, escape, bottom_escape, net, inelastic, ionization, excitation,
                        heating, subfloor, residual, input > 0 ? residual / input : NaN,
                        input > 0 ? escape / input : NaN, ionpairs,
                        inelastic_by_species, interval)
end

# Budget of one flux snapshot `Ie`, shape [n_z, n_μ, n_E].
function energy_budget_snapshot(model, Ie; verbose::Bool = true,
                                subfloor_factors = budget_subfloor_factors(model))
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
    axes(Ie) == (eachindex(z), eachindex(μ), eachindex(E)) ||
        throw(ArgumentError("Ie snapshot $(size(Ie)) does not match the model grid " *
                            "(n_z, n_μ, n_E) = ($n_z, $n_μ, $n_E)"))
    length(subfloor_factors.collisions) == length(model.species) &&
        all(size(f) == size(sp.cross_sections)
            for (f, sp) in zip(subfloor_factors.collisions, model.species)) &&
        size(subfloor_factors.coulomb) == (n_z, n_E) ||
        throw(DimensionMismatch("subfloor_factors do not match the model; compute them " *
                                "with budget_subfloor_factors(model)"))

    # ---- Boundary energy fluxes ----------------------------------------------------------
    # Vertical energy flux = Σ_beam Σ_E Ie·E·|μ| (matches field_aligned_beam_norm), taken at
    # the two points the transport operator's column sum telescopes to: the top row, and the
    # point above the absorbing bottom row, whose downward flux leaves the column.
    input = 0.0; escape = 0.0; bottom_escape = 0.0
    i_top = lastindex(Ie, 1)
    i_bottom = firstindex(Ie, 1) + 1
    for iμ in axes(Ie, 2), iE in axes(Ie, 3)
        μ[iμ] < 0 || continue
        bottom_escape += Ie[i_bottom, iμ, iE] * E[iE] * abs(μ[iμ])
    end
    for iμ in axes(Ie, 2), iE in axes(Ie, 3)
        fe = Ie[i_top, iμ, iE] * E[iE] * abs(μ[iμ])
        if μ[iμ] < 0
            input += fe
        else
            escape += fe
        end
    end

    # ---- Weighted omnidirectional flux for the volumetric (deposition / heating) terms ---
    # Every volumetric term is linear in the local flux, so each beam's quadrature weight is
    # folded into `Ie_omni` once. The profiles below are therefore already weighted, and their
    # plain sums over altitude are the column integrals.
    w = column_weights(s)
    Ie_omni = zeros(n_z, n_E)                              # [n_z, n_E], weight-folded
    for iμ in axes(Ie, 2)
        w_μ = μ[iμ] < 0 ? w.down : w.up
        for iE in axes(Ie, 3), iz in axes(Ie, 1)
            Ie_omni[iz, iE] += w_μ[iz] * Ie[iz, iμ, iE]
        end
    end

    # ---- Inelastic energy deposition (∫ along the field line) ----------------------------
    # Σ_s w_s Σ_sp n Σ_levels threshold·Σ_E Ie_omni·σ. The sub-floor term is the same sum with
    # the threshold replaced by `F`, the energy one collision takes below the grid.
    dep_profile_total = zeros(n_z)      # weighted energy deposition [eV m⁻² s⁻¹]
    ion_profile       = zeros(n_z)      # → ionization
    exc_profile       = zeros(n_z)      # → excitation
    ionpair_profile   = zeros(n_z)      # weighted ion-pair production [m⁻² s⁻¹]
    inelastic_by_species = Pair{String,Float64}[]
    subfloor = 0.0
    for (sp, F) in zip(model.species, subfloor_factors.collisions)
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
                acc_subfloor = 0.0
                for iE in eachindex(E)
                    acc += Ie_omni[iz, iE] * σ[lvl, iE]
                    acc_subfloor += Ie_omni[iz, iE] * σ[lvl, iE] * F[lvl, iE]
                end
                subfloor += dens[iz] * acc_subfloor
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

    # ---- Thermal-electron heating -------------------------------------------------------
    heating_profile = calculate_heating_rate(z, [0.0], reshape(Ie_omni, n_z, 1, n_E),
                                             E, ne, Te)[:, 1]
    heating = sum(heating_profile)

    # ---- Sub-floor energy of the Coulomb loss --------------------------------------------
    for iE in axes(Ie_omni, 2), iz in axes(Ie_omni, 1)
        subfloor += subfloor_factors.coulomb[iz, iE] * Ie_omni[iz, iE]
    end

    budget = close_budget(input, escape, bottom_escape, inelastic, ionization, excitation,
                          heating, subfloor, ionpairs, inelastic_by_species, nothing)
    verbose && (show(stdout, MIME"text/plain"(), budget); println())
    return budget
end

# Cell lengths for the column integral. The solvers use an upwind scheme, so at a given
# altitude a downward beam is differenced over the cell above it and an upward beam over the
# cell below it; each beam is integrated with the cell lengths its own equation uses. The
# first and last rows hold boundary conditions, not a balance, and get zero weight.
function column_weights(s)
    Base.require_one_based_indexing(s)
    Δs = abs.(diff(s))
    down = [0.0; Δs[2:end];     0.0]   # cell above, for μ < 0
    up   = [0.0; Δs[1:end - 1]; 0.0]   # cell below, for μ > 0
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

# ======================================================================================== #
#                 ENERGY OF THE ELECTRONS DEGRADED BELOW THE ENERGY GRID                    #
# ======================================================================================== #
#
# An inelastic collision or the Coulomb loss can take an electron below the lowest grid edge
# E_edges[1], out of the suprathermal population. The factors below give the energy those
# electrons carry, per collision or per unit Le·Ie, from the kinematics of each process
# alone.

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
  primary energies of the bin above `I`. The result is divided by the bin's ionization
  events (`event_count`), as the cascading matrices are. Double ionization uses the joint
  densities of the cascading calculation. The bin holding `I` has no cascading row and all
  of its products fall below the lowest edge, so its factor is the mean excess energy
  `(E_edges[iE+1] - I)/2` of the primaries above `I`.

These are the solver's own conventions for the two kinds of channel.
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
    # The law is evaluated millions of times, so it is unwrapped once and entered through a
    # single invokelatest.
    law = runtime_law(cache.spec.secondary_law)
    Base.invokelatest(fill_ionization_subfloor_energy!, factors, law, levels, ionizing_rows,
                      E_edges)

    # Per ionization event. The bin holding the threshold `I` has no cascading row:
    # `initialize!` guarantees that all of its products fall below the lowest edge, so each
    # event there carries its whole excess energy `E_p - I` below the floor, on average
    # half the bin's excess above `I`. A bin above the threshold with a cross section but
    # no events is rejected at `initialize!` as well.
    for l in ionizing_rows, iE in eachindex(ΔE)
        I = levels[l, 1]
        events = event_count(cache, iE, I)
        if events > 0
            factors[l, iE] /= events
        elseif σ[l, iE] > 0
            E_edges[iE] < I || throw(ArgumentError(
                "$(sp.name): ionizing channel at $(I) eV has cross section $(σ[l, iE]) m² \
                 in energy bin $(iE) but no ionization events in its cascading matrices"))
            factors[l, iE] = max(E_edges[iE + 1] - I, 0.0) / 2
        else
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
# the Coulomb factors times the loss function Le, [n_z, n_E].
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
