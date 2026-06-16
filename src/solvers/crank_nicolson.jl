# Optimized Crank-Nicolson scheme using direct nzval modification.
# This version avoids allocations by reusing sparse matrix structures for both
# Mlhs and Mrhs, writing physics values directly via pre-computed index arrays.

using KLU: klu, klu!
using LinearAlgebra: ldiv!, mul!
using SparseArrays: spdiagm

"""
    Crank_Nicolson!(Ie, t, model, v, matrices, iE, Ie_top, I0, cache)

Solve the time-dependent electron transport equation for energy level `iE`
using the Crank-Nicolson implicit scheme.

On the **first call** the sparse matrix structures, nzval index arrays, and
operator diagonals are computed and stored in `cache`.  On subsequent calls
only the numerical values in `Mlhs.nzval` / `Mrhs.nzval` are updated (zero
allocations on the hot path).

# Mathematical Background

The time-dependent electron transport equation is:
```
∂Ie/∂t + μ·v ∂Ie/∂z + A·Ie − ∫B·Ie′ dΩ′ − D·∇²Ie = Q
```

Crank-Nicolson gives second-order accuracy in time:
```
Mlhs · Ie^(n+1)  =  Mrhs · Ie^n  +  (Q^(n+1) + Q^n)/2
```

with
```
Mlhs = 1/(v·Δt) + μ·Ddz/2 + A/2 − B/2 − D·Ddiffusion
Mrhs = 1/(v·Δt) − μ·Ddz/2 − A/2 + B/2 + D·Ddiffusion
```

Both matrices share the same block structure as the steady-state system:
```
┌─────────┬─────────┬─────────┐
│ Block   │ Block   │ Block   │  Each block is n_z × n_z
│ (1,1)   │ (1,2)   │ (1,3)   │
├─────────┼─────────┼─────────┤
│ Block   │ Block   │ Block   │  Off-diagonal: angular scattering
│ (2,1)   │ (2,2)   │ (2,3)   │
├─────────┼─────────┼─────────┤
│ Block   │ Block   │ Block   │  Diagonal: transport + diffusion
│ (3,1)   │ (3,2)   │ (3,3)   │
└─────────┴─────────┴─────────┘
```

# Arguments
- `Ie`: pre-allocated output array [m⁻² s⁻¹], size `(n_z * n_angle, n_t)`
- `t`: time grid [s]
- `model`: `AuroraModel` (`s_field` and `pitch_angle_grid.μ_center` are used)
- `v`: electron velocity [km/s]
- `matrices::TransportMatrices`: container with `A`, `B`, `D`, `Q`, `Ddiffusion`
- `iE`: current energy index
- `Ie_top`: boundary condition at top [m⁻² s⁻¹] at each time step
- `I0`: initial condition [m⁻² s⁻¹]
- `cache`: `SolverCache` storing `Mlhs`, `Mrhs`, indices, `op_diags`, `KLU`
"""
function Crank_Nicolson!(Ie, t, model::AuroraModel, v, matrices, iE, Ie_top, I0, cache)
    z = model.s_field
    μ = model.pitch_angle_grid.μ_center
    n_z = length(z)
    n_angle = length(μ)

    # Extract physics data for this energy level
    A = matrices.A
    B = matrices.B
    D = @view(matrices.D[iE, :])
    Q_slice = @view(matrices.Q[:, :, iE])
    Ddiffusion = matrices.Ddiffusion

    # Temporal coefficient (scalar — all altitudes have the same 1/(v·Δt))
    dt = t[2] - t[1]
    ddt = 1.0 / (v * dt)

    # ── First call : build sparsity patterns, index maps, operator diags ──
    if !cache.initialized
        Ddz_Up, Ddz_Down = build_spatial_operators(z; half_weight = true)
        cache.Mlhs, cache.Mrhs = create_transport_sparsity_pattern(
            n_z, n_angle, μ, D, Ddiffusion; include_rhs = true)
        cache.indices_lhs = extract_nzval_indices(cache.Mlhs, n_z, n_angle)
        cache.indices_rhs = extract_nzval_indices(cache.Mrhs, n_z, n_angle)
        cache.op_diags    = extract_operator_diagonals(Ddz_Up, Ddz_Down, Ddiffusion)
        cache.rhs         = Vector{Float64}(undef, n_z * n_angle)
    end

    # ── Update matrix values (fast, no allocations) ──
    update_crank_nicolson_matrices!(cache.Mlhs, cache.Mrhs,
                                    cache.indices_lhs, cache.indices_rhs,
                                    A, B, D, ddt, cache.op_diags, μ, n_z)

    # ── Boundary indices ──
    index_bottom = 1:n_z:(n_angle * n_z)
    index_top    = n_z:n_z:(n_angle * n_z)

    # ── Initial condition ──
    Ie[:, 1] .= I0
    Ie[index_bottom, 1] .= 0.0
    Ie[index_top,    1] .= @view(Ie_top[:, 1])
    current = @view(Ie[:, 1])
    rhs = similar(current)

    # ── Factorise / re-factorise ──
    if !cache.initialized
        cache.KLU = klu(cache.Mlhs)
        cache.initialized = true
    else
        klu!(cache.KLU, cache.Mlhs)
    end

    # ── Time-stepping loop ──
    for i_t in 1:(length(t) - 1)
        next = @view Ie[:, i_t + 1]

        # Crank-Nicolson step:  Mlhs · Ie^(n+1)  =  Mrhs · Ie^n  +  Q
        mul!(rhs, cache.Mrhs, current)
        @views @. next = rhs + 0.5 * (Q_slice[:, i_t] + Q_slice[:, i_t + 1])
        next[index_bottom] .= 0.0                           # bottom BC
        next[index_top]    .= @view(Ie_top[:, i_t + 1])     # top BC
        ldiv!(cache.KLU, next)

        current = next
    end

    # Check for negative values (if it happens we have a problem) and clamp to zero
    if any(Ie .< 0)
        worst = minimum(Ie)
        worst_index = argmin(Ie)
        angle_index = (worst_index[1] - 1) ÷ n_z + 1
        z_index = mod1(worst_index[1], n_z)
        time_index = worst_index[2]
        energy = model.energy_grid.E_centers[iE]
        height = z[z_index]
        pitch_angle = acosd(μ[angle_index])

        @debug "Negative fluxes detected and clamped to zero ($(count(Ie .< 0)) values, worst = $(worst), " *
        "at [height = $(height) m, time = $(t[time_index]) s, energy = $(energy) eV, pitch_angle = $(pitch_angle) deg])"
        Ie[Ie .< 0] .= 0
    end

    return nothing
end


# ──────────────────────────────────────────────────────────────────────────────
# Matrix value update
# ──────────────────────────────────────────────────────────────────────────────

"""
    update_crank_nicolson_matrices!(Mlhs, Mrhs, idx_lhs, idx_rhs,
                                    A, B, D, ddt, op, μ, n_z)

Fill both `Mlhs` and `Mrhs` with the Crank-Nicolson operator values using the
pre-computed `BlockIndices` arrays and dense `OperatorDiagonals`.

`ddt` is the scalar `1/(v·Δt)` (constant for all altitudes).

The physics formulas (per stream direction) are:
```
Mlhs =  ddt  +  μ·Ddz  +  A/2  −  B/2  −  D·Ddiffusion
Mrhs =  ddt  −  μ·Ddz  −  A/2  +  B/2  +  D·Ddiffusion
```
where `Ddz` already contains the `/2` factor (built with `half_weight=true`).
"""
function update_crank_nicolson_matrices!(Mlhs, Mrhs, idx_lhs, idx_rhs,
                                         A, B, D, ddt::Float64,
                                         op::OperatorDiagonals, μ, n_z)
    n_angle = length(μ)
    nz_lhs = Mlhs.nzval
    nz_rhs = Mrhs.nzval
    interior = 2:(n_z - 1)

    for i1 in 1:n_angle
        for i2 in 1:n_angle
            bl = idx_lhs[i1, i2]
            br = idx_rhs[i1, i2]
            B_tmp = @view B[:, i1, i2]

            if i1 != i2
                # ── Off-diagonal blocks: scattering coupling ── #
                #   Mlhs: -B/2,   Mrhs: +B/2
                @views nz_lhs[bl.diag] .= .-B_tmp[interior] ./ 2
                @views nz_rhs[br.diag] .=   B_tmp[interior] ./ 2
            else
                # ── Diagonal blocks: transport + loss + diffusion ── #
                nz_lhs[bl.bc_first] = 1.0           # bottom boundary

                μ_i = μ[i1]
                D_i = D[i1]
                A_half = @view A[interior]           # used as A/2 below

                if μ_i < 0   # ── downward streams ──

                    # Main diagonal
                    #   Mlhs: ddt + μ·Ddz + A/2 - B/2 - D·Ddiff
                    #   Mrhs: ddt - μ·Ddz - A/2 + B/2 + D·Ddiff
                    @views nz_lhs[bl.diag] .= (ddt .+ μ_i .* op.Ddz_Down_diag[interior]
                                               .+ A_half ./ 2 .- B_tmp[interior] ./ 2
                                               .- D_i .* op.Ddiff_diag[interior])
                    @views nz_rhs[br.diag] .= (ddt .- μ_i .* op.Ddz_Down_diag[interior]
                                               .- A_half ./ 2 .+ B_tmp[interior] ./ 2
                                               .+ D_i .* op.Ddiff_diag[interior])

                    # Super-diagonal: μ·Ddz_super - D·Ddiff_super  /  negated
                    @views nz_lhs[bl.super] .= ( μ_i .* op.Ddz_Down_super[interior]
                                                .- D_i .* op.Ddiff_super[interior])
                    @views nz_rhs[br.super] .= (-μ_i .* op.Ddz_Down_super[interior]
                                                .+ D_i .* op.Ddiff_super[interior])

                    # Sub-diagonal (diffusion only, may be empty)
                    if !isempty(bl.sub)
                        @views nz_lhs[bl.sub] .= .- D_i .* op.Ddiff_sub[interior .- 1]
                        @views nz_rhs[br.sub] .=    D_i .* op.Ddiff_sub[interior .- 1]
                    end

                    # Top boundary
                    nz_lhs[bl.bc_last] = 1.0

                else         # ── upward streams ──

                    # Main diagonal
                    @views nz_lhs[bl.diag] .= (ddt .+ μ_i .* op.Ddz_Up_diag[interior]
                                               .+ A_half ./ 2 .- B_tmp[interior] ./ 2
                                               .- D_i .* op.Ddiff_diag[interior])
                    @views nz_rhs[br.diag] .= (ddt .- μ_i .* op.Ddz_Up_diag[interior]
                                               .- A_half ./ 2 .+ B_tmp[interior] ./ 2
                                               .+ D_i .* op.Ddiff_diag[interior])

                    # Sub-diagonal: μ·Ddz_sub - D·Ddiff_sub  /  negated
                    @views nz_lhs[bl.sub] .= ( μ_i .* op.Ddz_Up_sub[interior .- 1]
                                              .- D_i .* op.Ddiff_sub[interior .- 1])
                    @views nz_rhs[br.sub] .= (-μ_i .* op.Ddz_Up_sub[interior .- 1]
                                              .+ D_i .* op.Ddiff_sub[interior .- 1])

                    # Super-diagonal (diffusion only, may be empty)
                    if !isempty(bl.super)
                        @views nz_lhs[bl.super] .= .- D_i .* op.Ddiff_super[interior]
                        @views nz_rhs[br.super] .=    D_i .* op.Ddiff_super[interior]
                    end

                    # Top boundary: ∂Ie/∂z = 0  →  [-1, 1]
                    nz_lhs[bl.bc_last_sub] = -1.0
                    nz_lhs[bl.bc_last]     =  1.0
                end
            end
        end
    end

    return Mlhs, Mrhs
end
