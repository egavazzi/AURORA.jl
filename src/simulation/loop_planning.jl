function CFL_criteria(duration, dt, z, v, CFL_number=64)
    # The Courant-Freidrichs-Lewy (CFL) number normally has to be small (<4) to ensure numerical
    # stability. However, as a Crank-Nicolson scheme is always stable, we can take a bigger CFL. We
    # should be careful about numerical accuracy though.
    # For Gaussian inputs (or similar), it seems that the CFL can be set to 64 without major effects
    # on the results, while reducing computational time tremendously
    dz = z[2] - z[1]
    # Maximum dt that satisfies the CFL condition
    dt_max = CFL_number * dz / v
    # Refinement factor: how many internal steps per save interval.
    # dt_internal = dt / CFL_factor, which is exact by construction.
    CFL_factor = ceil(Int, dt / dt_max)
    # Number of save intervals (validated upstream to be an exact integer multiple of dt)
    n_save = round(Int, duration / dt)
    # Coarse save grid: the time points we want to write to disk.
    t_save = range(0.0, n_save * dt; length = n_save + 1)
    # Fine internal grid: CFL_factor sub-steps per save interval.
    t_internal = range(0.0, n_save * dt; length = n_save * CFL_factor + 1)

    return t_internal, t_save, CFL_factor
end


# ======================================================================================== #
#                       Memory estimation and automatic n_loop calculation                #
# ======================================================================================== #


"""
    calculate_n_loop(n_z, n_μ, n_E, n_t, N_neutrals, CFL_factor, n_save;
                     max_memory_gb=8, verbose=true) → Int

Choose the number of loops (`n_loop`) for a time-dependent simulation so that the peak
memory stays within `max_memory_gb`.

Using the split returned by [`estimate_simulation_memory`](@ref), the peak memory of an
`n_loop`-loop run is ≈ `fixed_gb + scaling_gb / n_loop`. The smallest `n_loop` satisfying
`fixed_gb + scaling_gb / n_loop ≤ max_memory_gb` is therefore

    n_loop = ceil(scaling_gb / (max_memory_gb - fixed_gb))

The result is capped at `n_save`, because a loop cannot hold fewer than one save interval.
If even one save interval per loop (or the fixed allocations alone) exceeds the budget, the
maximum split `n_loop = n_save` is returned together with a warning.

# Arguments
- `n_z`, `n_μ`, `n_E`: altitude, pitch-angle and energy grid sizes
- `n_t`: number of internal time steps over the *whole* simulation (full CFL-refined grid)
- `N_neutrals`: number of neutral species
- `CFL_factor`: internal sub-steps per save interval
- `n_save`: total number of save intervals

# Keyword Arguments
- `max_memory_gb`: memory budget (GB). Defaults to 8 GB.
- `verbose::Bool=true`: if true, print a short summary of the calculation

# Returns
- `n_loop::Int`: recommended number of loops
"""
function calculate_n_loop(n_z, n_μ, n_E, n_t, N_neutrals, CFL_factor, n_save;
                          max_memory_gb=8, verbose=true)
    mem = estimate_simulation_memory(n_z, n_μ, n_E, n_t, N_neutrals, CFL_factor)

    budget_unmet = false
    if mem.total_gb <= max_memory_gb
        n_loop = 1                                  # everything fits in a single loop
    elseif max_memory_gb <= mem.fixed_gb
        # The split-independent allocations already exceed the budget; splitting cannot
        # help, so fall back to the finest useful split.
        n_loop = n_save
        budget_unmet = true
    else
        n_loop = ceil(Int, mem.scaling_gb / (max_memory_gb - mem.fixed_gb))
        if n_loop > n_save                          # cannot split finer than 1 save/loop
            n_loop = n_save
            budget_unmet = true
        end
    end

    if verbose
        println("Automatic n_loop calculation:")
        println("  Grid dimensions:")
        println("    Altitude points (n_z):       $n_z")
        println("    Pitch-angle beams (n_μ):     $n_μ")
        println("    Time steps (n_t):            $n_t")
        println("    Energy points (n_E):         $n_E")
        println("  ─────────────────────────────────────────")
        println("  Fixed memory (not split):      $(round(mem.fixed_gb, digits=3)) GB")
        println("  Time-scaling memory (1 loop):  $(round(mem.scaling_gb, digits=3)) GB")
        println("  Single-loop total:             $(round(mem.total_gb, digits=3)) GB")
        println("  Total RAM detected:            $(round(Sys.total_memory() / 1e9, digits=2)) GB")
        println("  Max memory limit:              $(round(max_memory_gb, digits=2)) GB")
        println("  ─────────────────────────────────────────")
        println("  Calculated n_loop:             $n_loop")
        println()
    end

    if budget_unmet && verbose
        @warn "Cannot keep the peak memory within max_memory_gb = $(max_memory_gb) GB " *
              "even with n_loop = $n_save (one save interval per loop). " *
              "Estimated peak: $(round(mem.fixed_gb + mem.scaling_gb / n_save, digits=2)) GB."
    end

    return n_loop
end

"""
    check_n_loop(n_loop, n_z, n_μ, n_E, n_t, N_neutrals, CFL_factor; verbose=true)

Verify that running with `n_loop` loops fits in the detected RAM. Errors if the estimated
peak memory exceeds the available RAM (always), and warns if it exceeds half of it (only
when `verbose`).
"""
function check_n_loop(n_loop, n_z, n_μ, n_E, n_t, N_neutrals, CFL_factor; verbose::Bool=true)
    total_ram_gb = Sys.total_memory() / 1e9
    mem = estimate_simulation_memory(n_z, n_μ, n_E, n_t, N_neutrals, CFL_factor)
    peak_gb = mem.fixed_gb + mem.scaling_gb / n_loop

    if peak_gb > total_ram_gb
        # Smallest n_loop that would fit (if the fixed part alone already exceeds RAM,
        # no split can help).
        min_n_loop = mem.fixed_gb >= total_ram_gb ? nothing :
                     ceil(Int, mem.scaling_gb / (total_ram_gb - mem.fixed_gb))
        suggestion = isnothing(min_n_loop) ?
            "The fixed allocations alone ($(round(mem.fixed_gb, digits=2)) GB) exceed the " *
            "available RAM; reduce the grid size, duration or energy range." :
            "Increase `n_loop` to at least $min_n_loop, or decrease `max_memory_gb`."
        error("n_loop = $n_loop is too small to fit in the available RAM.\n" *
              "Estimated peak memory: $(round(peak_gb, digits=2)) GB\n" *
              "Detected RAM:          $(round(total_ram_gb, digits=2)) GB\n" *
              suggestion)
    elseif verbose && peak_gb > 0.5 * total_ram_gb
        @warn "n_loop = $n_loop may lead to high memory usage.\n" *
              "Estimated peak memory: $(round(peak_gb, digits=2)) GB\n" *
              "Detected RAM:          $(round(total_ram_gb, digits=2)) GB\n" *
              "Consider increasing `n_loop` or decreasing `max_memory_gb`."
    end
end

"""
    estimate_simulation_memory(n_z, n_μ, n_E, n_t, N_neutrals, CFL_factor)
        → (; fixed_gb, scaling_gb, total_gb)

Estimate the simulation's peak memory footprint (in decimal GB) by summing the actual
working arrays allocated in [`build_simulation_workspace`](@ref).

The footprint is split into two parts:

- `scaling_gb`: arrays laid out along the internal time grid — `workspace.Ie`, `matrices.Q`,
  the energy-degradation working arrays, and the subsampled `workspace.Ie_save` buffer.
  Evaluated here at the *full* internal length `n_t`; splitting the run into `n_loop`
  loops shrinks this contribution to ≈ `scaling_gb / n_loop`.
- `fixed_gb`: arrays whose size does **not** depend on the loop split — most importantly
  `workspace.Ie_top`, which always stores the entire input flux over the full time grid, plus
  `workspace.I0`, the transport/scattering matrices and the sparse solver factorisation.

`total_gb = fixed_gb + scaling_gb` is the single-loop (`n_loop = 1`) footprint, so the peak
memory of an `n_loop`-loop run is ≈ `fixed_gb + scaling_gb / n_loop`.
"""
function estimate_simulation_memory(n_z::Int, n_μ::Int, n_E::Int, n_t::Int,
                                    N_neutrals::Int, CFL_factor::Int)
    F = 8            # bytes per Float64
    GB = 1e9         # decimal GB, to match the user-facing `max_memory_gb`
    n_zμ = n_z * n_μ
    # Number of save intervals (the internal grid has n_t = n_save * CFL_factor + 1 points)
    n_save = (n_t - 1) ÷ CFL_factor

    # ── Arrays sized along the internal time grid (shrink by ≈ 1/n_loop when split) ──
    Ie      = n_zμ * n_t * n_E                  # workspace.Ie
    Q       = n_zμ * n_t * n_E                  # matrices.Q
    Ie_save = n_zμ * (n_save + 1) * n_E         # workspace.Ie_save (save cadence)
    # workspace.degradation working arrays, summed over the N_neutrals species:
    #   secondary_e_flux + primary_e_flux  → 2·N·(n_zμ·n_t)
    #   Ie_scatter                         → n_zμ·n_t
    #   ionization_source_sum              → n_z·n_t
    degradation = (2 * N_neutrals + 1) * n_zμ * n_t + n_z * n_t
    scaling_elems = Ie + Q + Ie_save + degradation

    # ── Arrays independent of the loop split ──
    Ie_top   = n_μ * n_t * n_E                  # workspace.Ie_top — full flux, never split
    I0       = n_zμ * n_E                       # workspace.I0
    # Transport matrices: A (n_z) + B (n_z·n_μ²)
    matrices = n_z + n_z * n_μ^2
    # Degradation energy spectra (2·N·n_E) + thermal-loss vector (n_z)
    misc     = 2 * N_neutrals * n_E + n_z
    # Sparse solver (Mlhs, Mrhs and the KLU factorisation). Both matrices share a block
    # structure with ≈ n_μ²·n_z stored values; KLU fill-in inflates this, so allow a
    # generous ×10 (still tiny next to Ie/Q because there is no n_E factor here).
    solver   = 10 * n_μ^2 * n_z
    fixed_elems = Ie_top + I0 + matrices + misc + solver

    scaling_gb = scaling_elems * F / GB
    fixed_gb   = fixed_elems   * F / GB
    return (; fixed_gb, scaling_gb, total_gb = fixed_gb + scaling_gb)
end
