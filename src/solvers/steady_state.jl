# Optimized steady state scheme using direct nzval modification.
# This version avoids allocations by reusing the sparse matrix structure and
# writing physics values directly into `nzval` via pre-computed index arrays.

using KLU: klu, klu!
using SparseArrays: spdiagm

# ──────────────────────────────────────────────────────────────────────────────
# Matrix value update
# ──────────────────────────────────────────────────────────────────────────────

"""
    update_steady_state_matrix!(Mlhs, indices, A, B, D, op, μ, n_z)

Fill the sparse matrix `Mlhs` with the steady-state transport operator values
using the pre-computed `indices` (a `Matrix{BlockIndices}`) and dense operator
diagonals `op` (`OperatorDiagonals`).

The physics formula for the system matrix is (per stream direction):
```
Mlhs = μ * Ddz  +  diag(A)  −  B  −  D * Ddiffusion
```
where `Ddz = Ddz_Down` for downward (μ < 0) or `Ddz_Up` for upward (μ > 0).
"""
function update_steady_state_matrix!(Mlhs, indices, A, B, D, op::OperatorDiagonals, μ, n_z)
    n_angle = length(μ)
    nzval = Mlhs.nzval
    interior = 2:(n_z - 1)

    for i1 in 1:n_angle
        for i2 in 1:n_angle
            bi = indices[i1, i2]
            B_tmp = @view B[:, i1, i2]

            if i1 != i2
                # ── Off-diagonal blocks: scattering coupling ──
                #   -B[interior]
                @views nzval[bi.diag] .= .-B_tmp[interior]
            else
                # ── Diagonal blocks: transport + loss + diffusion ──
                nzval[bi.bc_first] = 1.0          # bottom boundary

                μ_i = μ[i1]
                D_i = D[i1]

                if μ_i < 0   # ── downward streams ──
                    # Main diagonal:  μ*Ddz_Down + A - B - D*Ddiffusion
                    @views nzval[bi.diag] .= (μ_i .* op.Ddz_Down_diag[interior]
                                              .+ A[interior] .- B_tmp[interior]
                                              .- D_i .* op.Ddiff_diag[interior])

                    # Super-diagonal: μ*Ddz_Down_super - D*Ddiff_super
                    @views nzval[bi.super] .= (μ_i .* op.Ddz_Down_super[interior]
                                               .- D_i .* op.Ddiff_super[interior])

                    # Sub-diagonal (diffusion only, may be empty):  -D*Ddiff_sub
                    if !isempty(bi.sub)
                        @views nzval[bi.sub] .= .- D_i .* op.Ddiff_sub[interior .- 1]
                    end

                    # Top boundary: Ie = 0
                    nzval[bi.bc_last] = 1.0

                else         # ── upward streams ──
                    # Main diagonal:  μ*Ddz_Up + A - B - D*Ddiffusion
                    @views nzval[bi.diag] .= (μ_i .* op.Ddz_Up_diag[interior]
                                              .+ A[interior] .- B_tmp[interior]
                                              .- D_i .* op.Ddiff_diag[interior])

                    # Sub-diagonal: μ*Ddz_Up_sub - D*Ddiff_sub
                    @views nzval[bi.sub] .= (μ_i .* op.Ddz_Up_sub[interior .- 1]
                                             .- D_i .* op.Ddiff_sub[interior .- 1])

                    # Super-diagonal (diffusion only, may be empty):  -D*Ddiff_super
                    if !isempty(bi.super)
                        @views nzval[bi.super] .= .- D_i .* op.Ddiff_super[interior]
                    end

                    # Top boundary: ∂Ie/∂z = 0  →  [-1, 1]
                    nzval[bi.bc_last_sub] = -1.0
                    nzval[bi.bc_last]     =  1.0
                end
            end
        end
    end

    return Mlhs
end

# ──────────────────────────────────────────────────────────────────────────────
# Entry point
# ──────────────────────────────────────────────────────────────────────────────

"""
    steady_state_scheme!(Ie, model, matrices, iE, Ie_top, cache; first_iteration = false)

Solve the steady-state electron transport equation for energy level `iE`.

On the **first call** (`first_iteration = true`) the sparse matrix structure,
nzval index arrays, and operator diagonals are computed and stored in `cache`.
On subsequent calls only the numerical values in `Mlhs.nzval` are updated
(zero allocations on the hot path).

# Mathematical Background

The steady-state electron transport equation is:
```
μ ∂Ie/∂z + A·Ie − ∫B·Ie′ dΩ′ = Q
```

After spatial discretization this becomes the linear system  `Mlhs · Ie = Q`  with:
```
Mlhs = μ·Ddz + diag(A) − B − D·Ddiffusion
```

The matrix has a block structure indexed by pitch-angle pairs `(i1, i2)`:
```
┌─────────┬─────────┬─────────┐
│ Block   │ Block   │ Block   │  Each block is n_z × n_z
│ (1,1)   │ (1,2)   │ (1,3)   │  (n_z = number of altitudes)
├─────────┼─────────┼─────────┤
│ Block   │ Block   │ Block   │  Off-diagonal (i1≠i2): scattering (−B)
│ (2,1)   │ (2,2)   │ (2,3)   │
├─────────┼─────────┼─────────┤
│ Block   │ Block   │ Block   │  Diagonal (i1=i2): transport + loss + diffusion
│ (3,1)   │ (3,2)   │ (3,3)   │
└─────────┴─────────┴─────────┘
```

# Arguments
- `Ie`: pre-allocated output array [m⁻² s⁻¹], size `n_z * n_angle`
- `model`: `AuroraModel` (`s_field` and `pitch_angle_grid.μ_center` are used)
- `matrices::TransportMatrices`: container with `A`, `B`, `D`, `Q`, `Ddiffusion`
- `iE`: current energy index
- `Ie_top`: boundary condition at top [m⁻² s⁻¹]
- `cache`: `SolverCache` storing `Mlhs`, `indices_lhs`, `op_diags`, `KLU`
- `first_iteration`: whether this is the first call (creates structure) or subsequent (reuses)
"""
function steady_state_scheme!(Ie, model::AuroraModel, matrices, iE, Ie_top, cache; first_iteration = false)
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

    # ── First iteration: build sparsity pattern, index map, operator diagonals ──
    if first_iteration
        Ddz_Up, Ddz_Down = build_spatial_operators(z)
        cache.Mlhs       = create_transport_sparsity_pattern(n_z, n_angle, μ, D, Ddiffusion)
        cache.indices_lhs = extract_nzval_indices(cache.Mlhs, n_z, n_angle)
        cache.op_diags    = extract_operator_diagonals(Ddz_Up, Ddz_Down, Ddiffusion)
    end

    # ── Update matrix values (fast, no allocations) ──
    update_steady_state_matrix!(cache.Mlhs, cache.indices_lhs, A, B, D,
                                cache.op_diags, μ, n_z)

    # ── Factorise / re-factorise ──
    if first_iteration
        cache.KLU = klu(cache.Mlhs)
    else
        klu!(cache.KLU, cache.Mlhs)
    end

    # ── Boundary indices ──
    index_bottom = 1:n_z:(n_angle * n_z)
    index_top    = n_z:n_z:(n_angle * n_z)

    # ── Set boundary conditions in the RHS vector ──
    Q_local = copy(Q_slice)
    Q_local[index_bottom] .= 0.0
    Q_local[index_top]    .= Ie_top

    # ── Solve  Mlhs · Ie = Q ──
    Ie .= cache.KLU \ Q_local

    # Check for negative values (if it happens we have a problem) and clamp to zero
    if any(Ie .< 0)
        @warn "Negative fluxes detected and clamped to zero ($(count(Ie .< 0)) values)"
        Ie[Ie .< 0] .= 0
    end

    return nothing
end
