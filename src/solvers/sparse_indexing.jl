# Shared sparse matrix indexing infrastructure for the transport solvers.
#
# The steady-state and Crank-Nicolson solvers both build sparse matrices with the
# same block structure (one n_z × n_z block per angle pair).  This module provides:
#
#   1. `BlockIndices`       – typed struct replacing Dict{Tuple{Symbol,Int},Int}
#   2. `OperatorDiagonals`  – dense vectors extracted from the sparse FD operators
#   3. `create_transport_sparsity_pattern`  – builds Mlhs (and optionally Mrhs)
#   4. `extract_nzval_indices`              – maps block positions → nzval indices
#   5. `extract_operator_diagonals`         – extracts dense diagonals from sparse ops

using LinearAlgebra: diag
using SparseArrays: SparseMatrixCSC, sparse, spdiagm

# ──────────────────────────────────────────────────────────────────────────────
# Index storage (replaces Dict-based mapping)
# ──────────────────────────────────────────────────────────────────────────────

"""
    BlockIndices

Pre-computed `nzval` indices for a single (i1, i2) block of the transport matrix.

For **off-diagonal blocks** (i1 ≠ i2) only `diag` is populated (scattering coupling),
and the boundary fields are zero.

For **diagonal blocks** (i1 = i2) the tridiagonal entries and boundary condition
indices are stored.  `bc_last_sub` is zero for downward streams (no sub-diagonal in
the last-row boundary condition).
"""
struct BlockIndices
    diag::Vector{Int}    # main diagonal,  interior points 2:(n_z-1)
    super::Vector{Int}   # super-diagonal, interior points (may be empty)
    sub::Vector{Int}     # sub-diagonal,   interior points (may be empty)
    bc_first::Int        # first-row boundary (0 for off-diagonal blocks)
    bc_last::Int         # last-row  boundary (0 for off-diagonal blocks)
    bc_last_sub::Int     # last-row  sub-diag (0 when absent)
end

# ──────────────────────────────────────────────────────────────────────────────
# Cached dense diagonals of the finite-difference operators
# ──────────────────────────────────────────────────────────────────────────────

"""
    OperatorDiagonals

Dense vectors extracted once from the sparse finite-difference operators `Ddz_Up`,
`Ddz_Down` and `Ddiffusion`.  Storing them avoids repeated CSC look-ups during the
value-update phase.

Fields with suffix `_diag` hold the main diagonal, `_sub` the sub-diagonal
(offset −1), and `_super` the super-diagonal (offset +1).
"""
struct OperatorDiagonals
    # Ddz_Up  (used for upward streams, μ > 0)
    Ddz_Up_diag::Vector{Float64}    # length n_z
    Ddz_Up_sub::Vector{Float64}     # length n_z−1

    # Ddz_Down (used for downward streams, μ < 0)
    Ddz_Down_diag::Vector{Float64}  # length n_z
    Ddz_Down_super::Vector{Float64} # length n_z−1

    # Ddiffusion (pitch-angle diffusion ∂²/∂z²)
    Ddiff_diag::Vector{Float64}     # length n_z
    Ddiff_sub::Vector{Float64}      # length n_z−1
    Ddiff_super::Vector{Float64}    # length n_z−1
end

"""
    extract_operator_diagonals(Ddz_Up, Ddz_Down, Ddiffusion)

Extract the tridiagonal entries of the three finite-difference operators into
dense vectors for fast element-wise access.
"""
function extract_operator_diagonals(Ddz_Up, Ddz_Down, Ddiffusion)
    return OperatorDiagonals(
        diag(Ddz_Up, 0),  diag(Ddz_Up, -1),
        diag(Ddz_Down, 0), diag(Ddz_Down, 1),
        diag(Ddiffusion, 0), diag(Ddiffusion, -1), diag(Ddiffusion, 1),
    )
end

# ──────────────────────────────────────────────────────────────────────────────
# Sparsity pattern construction
# ──────────────────────────────────────────────────────────────────────────────

"""
    create_transport_sparsity_pattern(n_z, n_angle, μ, D, Ddiffusion; include_rhs=false)

Build the sparse matrix with the correct sparsity structure for the transport
system.  The returned matrix contains placeholder values; the actual physics
values are filled in later by the update functions.

When `include_rhs=true` a second matrix `Mrhs` (same structure but without
boundary rows) is also returned – used by Crank-Nicolson.

Returns `Mlhs` or `(Mlhs, Mrhs)`.
"""
function create_transport_sparsity_pattern(n_z, n_angle, μ, D, Ddiffusion; include_rhs::Bool = false)
    max_nnz = n_angle * n_angle * 3 * n_z
    row_l, col_l, val_l = _alloc_coo(max_nnz)
    row_r, col_r, val_r = include_rhs ? _alloc_coo(max_nnz) : (Int[], Int[], Float64[])

    for i1 in 1:n_angle
        for i2 in 1:n_angle
            if i1 != i2
                _add_offdiagonal_block!(row_l, col_l, val_l, i1, i2, n_z)
                include_rhs && _add_offdiagonal_block!(row_r, col_r, val_r, i1, i2, n_z)
            else
                _add_diagonal_block_lhs!(row_l, col_l, val_l, i1, n_z, μ, D, Ddiffusion)
                include_rhs && _add_diagonal_block_rhs!(row_r, col_r, val_r, i1, n_z, μ, D, Ddiffusion)
            end
        end
    end

    N = n_z * n_angle
    Mlhs = sparse(row_l, col_l, val_l, N, N)
    if include_rhs
        Mrhs = sparse(row_r, col_r, val_r, N, N)
        return Mlhs, Mrhs
    end
    return Mlhs
end

function _alloc_coo(max_nnz)
    rows = Vector{Int}();  sizehint!(rows, max_nnz)
    cols = Vector{Int}();  sizehint!(cols, max_nnz)
    vals = Vector{Float64}(); sizehint!(vals, max_nnz)
    return rows, cols, vals
end

function _add_offdiagonal_block!(rows, cols, vals, i1, i2, n_z)
    offset_row = (i1 - 1) * n_z
    offset_col = (i2 - 1) * n_z
    for i in 2:(n_z - 1)
        push!(rows, offset_row + i)
        push!(cols, offset_col + i)
        push!(vals, 0.0)
    end
end

function _add_diagonal_block_lhs!(rows, cols, vals, i1, n_z, μ, D, Ddiffusion)
    offset = (i1 - 1) * n_z

    # First row boundary condition
    push!(rows, offset + 1); push!(cols, offset + 1); push!(vals, 1.0)

    if μ[i1] < 0    # downward fluxes
        for i in 2:(n_z - 1)
            push!(rows, offset + i); push!(cols, offset + i);     push!(vals, 0.0)  # diag
            push!(rows, offset + i); push!(cols, offset + i + 1); push!(vals, 0.0)  # super
            if D[i1] != 0 && Ddiffusion[i, i - 1] != 0                              # sub (diffusion)
                push!(rows, offset + i); push!(cols, offset + i - 1); push!(vals, 0.0)
            end
        end
        # Last row boundary
        push!(rows, offset + n_z); push!(cols, offset + n_z); push!(vals, 1.0)
    else             # upward fluxes
        for i in 2:(n_z - 1)
            push!(rows, offset + i); push!(cols, offset + i);     push!(vals, 0.0)  # diag
            push!(rows, offset + i); push!(cols, offset + i - 1); push!(vals, 0.0)  # sub
            if D[i1] != 0 && Ddiffusion[i, i + 1] != 0                              # super (diffusion)
                push!(rows, offset + i); push!(cols, offset + i + 1); push!(vals, 0.0)
            end
        end
        # Last row boundary
        push!(rows, offset + n_z); push!(cols, offset + n_z - 1); push!(vals, -1.0)
        push!(rows, offset + n_z); push!(cols, offset + n_z);     push!(vals, 1.0)
    end
end

function _add_diagonal_block_rhs!(rows, cols, vals, i1, n_z, μ, D, Ddiffusion)
    offset = (i1 - 1) * n_z
    # No boundary rows in Mrhs

    if μ[i1] < 0    # downward fluxes
        for i in 2:(n_z - 1)
            push!(rows, offset + i); push!(cols, offset + i);     push!(vals, 0.0)  # diag
            push!(rows, offset + i); push!(cols, offset + i + 1); push!(vals, 0.0)  # super
            if D[i1] != 0 && Ddiffusion[i, i - 1] != 0                              # sub (diffusion)
                push!(rows, offset + i); push!(cols, offset + i - 1); push!(vals, 0.0)
            end
        end
    else             # upward fluxes
        for i in 2:(n_z - 1)
            push!(rows, offset + i); push!(cols, offset + i);     push!(vals, 0.0)  # diag
            push!(rows, offset + i); push!(cols, offset + i - 1); push!(vals, 0.0)  # sub
            if D[i1] != 0 && Ddiffusion[i, i + 1] != 0                              # super (diffusion)
                push!(rows, offset + i); push!(cols, offset + i + 1); push!(vals, 0.0)
            end
        end
    end
end

# ──────────────────────────────────────────────────────────────────────────────
# nzval index extraction
# ──────────────────────────────────────────────────────────────────────────────

"""
    extract_nzval_indices(M, n_z, n_angle)

Walk the CSC structure of sparse matrix `M` and return a `Matrix{BlockIndices}`
(size `n_angle × n_angle`) that maps each block's tridiagonal positions to their
index in `M.nzval`.

This is computed once after sparsity-pattern creation, then used on every energy
step to write values directly into `nzval` without any Dict look-up.
"""
function extract_nzval_indices(M::SparseMatrixCSC, n_z::Int, n_angle::Int)
    # Build reverse lookup: (row, col) → nzval index
    nz_lookup = Dict{Tuple{Int,Int}, Int}()
    sizehint!(nz_lookup, length(M.nzval))
    colptr = M.colptr
    rowval = M.rowval
    for col in 1:size(M, 2)
        for idx in colptr[col]:(colptr[col + 1] - 1)
            nz_lookup[(rowval[idx], col)] = idx
        end
    end

    indices = Matrix{BlockIndices}(undef, n_angle, n_angle)

    for i1 in 1:n_angle
        for i2 in 1:n_angle
            offset_row = (i1 - 1) * n_z
            offset_col = (i2 - 1) * n_z

            if i1 != i2
                # Off-diagonal: only main diagonal entries
                d = _collect_indices(nz_lookup, offset_row, offset_col, 0, n_z)
                indices[i1, i2] = BlockIndices(d, Int[], Int[], 0, 0, 0)
            else
                offset = (i1 - 1) * n_z
                d     = _collect_indices(nz_lookup, offset, offset,  0, n_z)
                sup   = _collect_indices(nz_lookup, offset, offset, +1, n_z)
                sub   = _collect_indices(nz_lookup, offset, offset, -1, n_z)

                bc_first   = get(nz_lookup, (offset + 1,    offset + 1),        0)
                bc_last    = get(nz_lookup, (offset + n_z,  offset + n_z),      0)
                bc_last_sub = get(nz_lookup, (offset + n_z, offset + n_z - 1),  0)

                indices[i1, i2] = BlockIndices(d, sup, sub, bc_first, bc_last, bc_last_sub)
            end
        end
    end

    return indices
end

"""
Collect nzval indices for a given diagonal offset within the interior rows 2:(n_z-1).
`offset` is 0 for main diagonal, +1 for super, -1 for sub.
"""
function _collect_indices(nz_lookup, offset_row, offset_col, diag_offset, n_z)
    out = Int[]
    sizehint!(out, n_z - 2)
    for i in 2:(n_z - 1)
        idx = get(nz_lookup, (offset_row + i, offset_col + i + diag_offset), 0)
        if idx != 0
            push!(out, idx)
        end
    end
    return out
end

# ──────────────────────────────────────────────────────────────────────────────
# Spatial differentiation matrix helpers
# ──────────────────────────────────────────────────────────────────────────────

"""
    build_spatial_operators(z; half_weight=false)

Construct the upwind/downwind spatial differentiation matrices `Ddz_Up` and
`Ddz_Down` from the altitude grid `z`.

When `half_weight=true` (used by Crank-Nicolson), the coefficients are halved
so that the matrices represent `Ddz/2` directly.
"""
function build_spatial_operators(z; half_weight::Bool = false)
    h4diffu = [z[1] - (z[2] - z[1]); z]
    h4diffd = [z; z[end] + (z[end] - z[end - 1])]
    scale = half_weight ? 2.0 : 1.0
    Ddz_Up   = spdiagm(-1 => -1 ./ (scale .* diff(h4diffu[2:end])),
                         0 =>  1 ./ (scale .* diff(h4diffu[1:end])))
    Ddz_Down = spdiagm( 0 => -1 ./ (scale .* diff(h4diffd[1:end])),
                         1 =>  1 ./ (scale .* diff(h4diffd[1:(end - 1)])))
    return Ddz_Up, Ddz_Down
end
