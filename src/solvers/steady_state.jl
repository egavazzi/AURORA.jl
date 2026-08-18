# Optimized steady state scheme using direct nzval modification.
# This version avoids allocations by reusing the sparse matrix structure and
# writing physics values directly into `nzval` via pre-computed index arrays.

using KLU: klu, klu!
using LinearAlgebra: ldiv!
using SparseArrays: spdiagm

# ──────────────────────────────────────────────────────────────────────────────
# Matrix value update
# ──────────────────────────────────────────────────────────────────────────────

"""
    update_steady_state_matrix!(Mlhs, indices, A, B, op, μ, n_z)

Fill the sparse matrix `Mlhs` with the steady-state transport operator values
using the pre-computed `indices` (a `Matrix{BlockIndices}`) and dense operator
diagonals `op` (`OperatorDiagonals`).

The physics formula for the system matrix is (per stream direction):
```
Mlhs = μ * Ddz  +  diag(A)  −  B
```
where `Ddz = Ddz_Down` for downward (μ < 0) or `Ddz_Up` for upward (μ > 0).
"""
function update_steady_state_matrix!(Mlhs, indices, A, B, op::OperatorDiagonals, μ, n_z)
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
                # ── Diagonal blocks: transport + loss ──
                nzval[bi.bc_first] = 1.0          # bottom boundary

                μ_i = μ[i1]

                if μ_i < 0   # ── downward streams ──
                    # Main diagonal:  μ*Ddz_Down + A - B
                    @views nzval[bi.diag] .= (μ_i .* op.Ddz_Down_diag[interior]
                                              .+ A[interior] .- B_tmp[interior])

                    # Super-diagonal: μ*Ddz_Down_super
                    @views nzval[bi.super] .= μ_i .* op.Ddz_Down_super[interior]

                    # Top boundary: Ie = 0
                    nzval[bi.bc_last] = 1.0

                else         # ── upward streams ──
                    # Main diagonal:  μ*Ddz_Up + A - B
                    @views nzval[bi.diag] .= (μ_i .* op.Ddz_Up_diag[interior]
                                              .+ A[interior] .- B_tmp[interior])

                    # Sub-diagonal: μ*Ddz_Up_sub
                    @views nzval[bi.sub] .= μ_i .* op.Ddz_Up_sub[interior .- 1]

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
    steady_state_scheme!(Ie, model, matrices, iE, Ie_top, workspace)

Solve the steady-state electron transport equation for energy level `iE`.

On the **first call** the sparse matrix structure, nzval index arrays, and operator diagonals
are computed and stored in `workspace`. On subsequent calls only the numerical values in
`Mlhs.nzval` are updated (zero allocations on the hot path).

# Mathematical Background

The steady-state electron transport equation is:
```
μ ∂Ie/∂z + A·Ie − ∫B·Ie′ dΩ′ = Q
```

After spatial discretization this becomes the linear system  `Mlhs · Ie = Q`  with:
```
Mlhs = μ·Ddz + diag(A) − B
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
│ Block   │ Block   │ Block   │  Diagonal (i1=i2): transport + loss
│ (3,1)   │ (3,2)   │ (3,3)   │
└─────────┴─────────┴─────────┘
```

# Arguments
- `Ie`: pre-allocated output array [m⁻² s⁻¹], size `n_z * n_angle`
- `model`: `AuroraModel` (`s_field` and `pitch_angle_grid.μ_center` are used)
- `matrices::TransportMatrices`: container with `A`, `B`, `Q`
- `iE`: current energy index
- `Ie_top`: boundary condition at top [m⁻² s⁻¹]
- `workspace`: `SolverWorkspace` storing `Mlhs`, `indices_lhs`, `op_diags`, `KLU`
"""
function steady_state_scheme!(Ie, model::AuroraModel, matrices, iE, Ie_top, workspace)
    z = model.s_field
    μ = model.pitch_angle_grid.μ_center
    n_z = length(z)
    n_angle = length(μ)

    # Extract physics data for this energy level
    A = matrices.A
    B = matrices.B
    Q_slice = @view(matrices.Q[:, :, iE])

    # ── First call: build sparsity pattern, index map, operator diagonals ──
    if !workspace.initialized
        Ddz_Up, Ddz_Down = build_spatial_operators(z)
        workspace.Mlhs = create_transport_sparsity_pattern(n_z, n_angle, μ)
        workspace.indices_lhs = extract_nzval_indices(workspace.Mlhs, n_z, n_angle)
        workspace.op_diags = extract_operator_diagonals(Ddz_Up, Ddz_Down)
    end

    # ── Update matrix values (fast, no allocations) ──
    update_steady_state_matrix!(workspace.Mlhs, workspace.indices_lhs, A, B,
                                workspace.op_diags, μ, n_z)

    # ── Factorise / re-factorise ──
    if !workspace.initialized
        workspace.KLU = klu(workspace.Mlhs)
        workspace.initialized = true
    else
        klu!(workspace.KLU, workspace.Mlhs)
    end

    # ── Boundary indices ──
    index_bottom = 1:n_z:(n_angle * n_z)
    index_top    = n_z:n_z:(n_angle * n_z)

    # ── Set up the RHS in `Ie`, apply boundary conditions, and solve in place ──
    copyto!(Ie, Q_slice)
    Ie[index_bottom] .= 0.0
    Ie[index_top]    .= Ie_top
    ldiv!(workspace.KLU, Ie)

    # Check for negative values (if it happens we have a problem) and clamp to zero
    if any(<(0), Ie)
        @warn "Negative fluxes detected and clamped to zero ($(count(<(0), Ie)) values)"
        clamp!(Ie, 0.0, Inf)
    end

    return nothing
end
