# Optimized steady state scheme using direct nzval modification
# This version avoids allocations by reusing sparse matrix structure

using KLU: klu, klu!
using LinearAlgebra: Diagonal
using SparseArrays: spdiagm, sparse, sparse!, dropzeros!, findnz

"""
    create_steady_state_sparsity_pattern(n_z, n_angle, μ, D, Ddiffusion)

Create the sparsity pattern (structure) of the steady-state LHS matrix once.
This can be reused by only modifying nzval values, avoiding allocations.

Returns the sparse matrix with the correct structure.
"""
function create_steady_state_sparsity_pattern(n_z, n_angle, μ, D, Ddiffusion)
    row_l = Vector{Int}()
    col_l = Vector{Int}()
    val_l = Vector{Float64}()

    # Pre-allocate with estimated size to avoid reallocations
    max_nnz = n_angle * n_angle * 3 * n_z
    sizehint!(row_l, max_nnz)
    sizehint!(col_l, max_nnz)
    sizehint!(val_l, max_nnz)

    for i1 in 1:n_angle
        for i2 in 1:n_angle
            if i1 != i2
                # Off-diagonal blocks: diagonal elements
                offset_row = (i1 - 1) * n_z
                offset_col = (i2 - 1) * n_z
                for i in 2:(n_z - 1)
                    push!(row_l, offset_row + i)
                    push!(col_l, offset_col + i)
                    push!(val_l, 0.0)  # placeholder
                end
            else
                offset = (i1 - 1) * n_z

                # First row boundary condition
                push!(row_l, offset + 1)
                push!(col_l, offset + 1)
                push!(val_l, 1.0)

                if μ[i1] < 0    # downward fluxes
                    # Main tridiagonal part (skip first and last)
                    for i in 2:(n_z - 1)
                        # Diagonal
                        push!(row_l, offset + i)
                        push!(col_l, offset + i)
                        push!(val_l, 0.0)

                        # Super-diagonal
                        push!(row_l, offset + i)
                        push!(col_l, offset + i + 1)
                        push!(val_l, 0.0)

                        # Sub-diagonal (from diffusion if present)
                        if D[i1] != 0 && Ddiffusion[i, i - 1] != 0
                            push!(row_l, offset + i)
                            push!(col_l, offset + i - 1)
                            push!(val_l, 0.0)
                        end
                    end

                    # Last row boundary condition
                    push!(row_l, offset + n_z)
                    push!(col_l, offset + n_z)
                    push!(val_l, 1.0)

                else  # upward fluxes
                    # Main tridiagonal part (skip first and last)
                    for i in 2:(n_z - 1)
                        # Diagonal
                        push!(row_l, offset + i)
                        push!(col_l, offset + i)
                        push!(val_l, 0.0)

                        # Sub-diagonal
                        push!(row_l, offset + i)
                        push!(col_l, offset + i - 1)
                        push!(val_l, 0.0)

                        # Super-diagonal (from diffusion if present)
                        if D[i1] != 0 && Ddiffusion[i, i + 1] != 0
                            push!(row_l, offset + i)
                            push!(col_l, offset + i + 1)
                            push!(val_l, 0.0)
                        end
                    end

                    # Last row boundary condition
                    push!(row_l, offset + n_z)
                    push!(col_l, offset + n_z - 1)
                    push!(val_l, -1.0)

                    push!(row_l, offset + n_z)
                    push!(col_l, offset + n_z)
                    push!(val_l, 1.0)
                end
            end
        end
    end

    return sparse(row_l, col_l, val_l, n_z * n_angle, n_z * n_angle)
end

"""
    create_steady_state_nzval_mapping(Mlhs, n_z, n_angle)

Create a mapping from matrix block positions to nzval indices for efficient in-place updates.

Julia uses the CSC (Compressed Sparse Column) format for sparse matrices, with three arrays:
- colptr: Column pointers (which rows are in each column)
- rowval: Row indices of non-zero values
- nzval: The actual non-zero values
This function computes a mapping that tells us which index in `nzval` corresponds to each
matrix element we want to update. With this mapping, we can directly modify `nzval` in
place, avoiding expensive matrix reconstruction.

Returns a structured mapping where `mapping[i1, i2][(:type, position)]` gives the nzval index
for the specified matrix element in block (i1, i2).
"""
function create_steady_state_nzval_mapping(Mlhs, n_z, n_angle)
    # For each (i1, i2) block, store the nzval indices
    # Structure: mapping[i1, i2] = Dict with keys like (:diag, i), (:super, i), (:sub, i)
    mapping = Matrix{Dict{Tuple{Symbol,Int},Int}}(undef, n_angle, n_angle)

    nzval = Mlhs.nzval
    rowval = Mlhs.rowval
    colptr = Mlhs.colptr

    # Build reverse lookup: (row, col) -> nzval_idx
    nz_lookup = Dict{Tuple{Int,Int}, Int}()
    for col in 1:size(Mlhs, 2)
        for idx in colptr[col]:(colptr[col+1]-1)
            row = rowval[idx]
            nz_lookup[(row, col)] = idx
        end
    end

    # Now create structured mapping for each block
    for i1 in 1:n_angle
        for i2 in 1:n_angle
            block_map = Dict{Tuple{Symbol,Int},Int}()
            offset_row = (i1 - 1) * n_z
            offset_col = (i2 - 1) * n_z

            if i1 != i2
                # Off-diagonal: only diagonal elements
                for i in 2:(n_z - 1)
                    idx = get(nz_lookup, (offset_row + i, offset_col + i), nothing)
                    if idx !== nothing
                        block_map[(:diag, i)] = idx
                    end
                end
            else
                # Diagonal blocks
                offset = (i1 - 1) * n_z

                # Boundary conditions
                idx = get(nz_lookup, (offset + 1, offset + 1), nothing)
                if idx !== nothing
                    block_map[(:bc_first, 1)] = idx
                end

                # Interior points
                for i in 2:(n_z - 1)
                    # Main diagonal
                    idx = get(nz_lookup, (offset + i, offset + i), nothing)
                    if idx !== nothing
                        block_map[(:diag, i)] = idx
                    end

                    # Super-diagonal
                    idx = get(nz_lookup, (offset + i, offset + i + 1), nothing)
                    if idx !== nothing
                        block_map[(:super, i)] = idx
                    end

                    # Sub-diagonal
                    idx = get(nz_lookup, (offset + i, offset + i - 1), nothing)
                    if idx !== nothing
                        block_map[(:sub, i)] = idx
                    end
                end

                # Last row boundary
                idx = get(nz_lookup, (offset + n_z, offset + n_z), nothing)
                if idx !== nothing
                    block_map[(:bc_last, n_z)] = idx
                end
                idx = get(nz_lookup, (offset + n_z, offset + n_z - 1), nothing)
                if idx !== nothing
                    block_map[(:bc_last_sub, n_z)] = idx
                end
            end

            mapping[i1, i2] = block_map
        end
    end

    return mapping
end

"""
    update_steady_state_matrix!(Mlhs, mapping, A, B, D, Ddiffusion, Ddz_Up, Ddz_Down, μ, h_atm)

Update the sparse matrix values using pre-computed mapping.
This avoids all allocations by directly modifying nzval.
"""
function update_steady_state_matrix!(Mlhs, mapping, A, B, D, Ddiffusion, Ddz_Up, Ddz_Down, μ, h_atm)
    n_z = length(h_atm)
    n_angle = length(μ)
    nzval = Mlhs.nzval

    for i1 in 1:n_angle
        for i2 in 1:n_angle
            B_tmp = @view B[:, i1, i2]
            block_map = mapping[i1, i2]

            if i1 != i2
                # Off-diagonal blocks
                for i in 2:(n_z - 1)
                    idx = get(block_map, (:diag, i), nothing)
                    if idx !== nothing
                        nzval[idx] = -B_tmp[i]
                    end
                end
            else
                # Boundary conditions
                idx = get(block_map, (:bc_first, 1), nothing)
                if idx !== nothing
                    nzval[idx] = 1.0
                end

                if μ[i1] < 0    # downward fluxes
                    for i in 2:(n_z - 1)
                        # Diagonal
                        val = μ[i1] * Ddz_Down[i, i] + A[i] - B_tmp[i]
                        if D[i1] != 0 && Ddiffusion[i, i] != 0
                            val -= D[i1] * Ddiffusion[i, i]
                        end
                        idx = get(block_map, (:diag, i), nothing)
                        if idx !== nothing
                            nzval[idx] = val
                        end

                        # Super-diagonal
                        val = μ[i1] * Ddz_Down[i, i + 1]
                        if D[i1] != 0 && Ddiffusion[i, i + 1] != 0
                            val -= D[i1] * Ddiffusion[i, i + 1]
                        end
                        idx = get(block_map, (:super, i), nothing)
                        if idx !== nothing
                            nzval[idx] = val
                        end

                        # Sub-diagonal
                        if D[i1] != 0 && Ddiffusion[i, i - 1] != 0
                            idx = get(block_map, (:sub, i), nothing)
                            if idx !== nothing
                                nzval[idx] = -D[i1] * Ddiffusion[i, i - 1]
                            end
                        end
                    end

                    # Last row boundary
                    idx = get(block_map, (:bc_last, n_z), nothing)
                    if idx !== nothing
                        nzval[idx] = 1.0
                    end

                else  # upward fluxes
                    for i in 2:(n_z - 1)
                        # Diagonal
                        val = μ[i1] * Ddz_Up[i, i] + A[i] - B_tmp[i]
                        if D[i1] != 0 && Ddiffusion[i, i] != 0
                            val -= D[i1] * Ddiffusion[i, i]
                        end
                        idx = get(block_map, (:diag, i), nothing)
                        if idx !== nothing
                            nzval[idx] = val
                        end

                        # Sub-diagonal
                        val = μ[i1] * Ddz_Up[i, i - 1]
                        if D[i1] != 0 && Ddiffusion[i, i - 1] != 0
                            val -= D[i1] * Ddiffusion[i, i - 1]
                        end
                        idx = get(block_map, (:sub, i), nothing)
                        if idx !== nothing
                            nzval[idx] = val
                        end

                        # Super-diagonal
                        if D[i1] != 0 && Ddiffusion[i, i + 1] != 0
                            idx = get(block_map, (:super, i), nothing)
                            if idx !== nothing
                                nzval[idx] = -D[i1] * Ddiffusion[i, i + 1]
                            end
                        end
                    end

                    # Last row boundary
                    idx = get(block_map, (:bc_last_sub, n_z), nothing)
                    if idx !== nothing
                        nzval[idx] = -1.0
                    end
                    idx = get(block_map, (:bc_last, n_z), nothing)
                    if idx !== nothing
                        nzval[idx] = 1.0
                    end
                end
            end
        end
    end

    return Mlhs
end

"""
    steady_state_scheme_optimized!(Ie, h_atm, μ, matrices, iE, Ie_top, cache; first_iteration = false)

Optimized steady-state scheme using direct nzval modification.
This is an in-place version that modifies `Ie` directly to avoid allocations.
This version avoids allocations by reusing the sparse matrix structure.

On first iteration, creates the sparsity pattern and mapping which are stored in cache.
On subsequent iterations, only updates the nzval array directly.

# Mathematical Background

The steady-state electron transport equation is:
```
μ ∂Ie/∂z + A*Ie - ∫B*Ie'dΩ' = Q
```

After spatial discretization, this becomes a linear system of coupled equations:
```
[μ*Ddz + A - B - D*Ddiffusion] * Ie = Q
         ↑
        Mlhs (the system matrix)
```

Where:
- `Ddz = Ddz_Up` or `Ddz_Down` (depending on sign of μ): spatial derivative operator
- `A`: electron loss rate matrix (diagonal)
- `B`: scattering operator matrix (couples different angles)
- `D`: diffusion coefficient (diagonal in angle space)
- `Ddiffusion`: second derivative operator for pitch-angle diffusion
- `Q`: source term (energy degradation and ionization)

The resulting sparse matrix `Mlhs` has a block structure:
```
┌─────────┬─────────┬─────────┐
│ Block   │ Block   │ Block   │  Each block is n_z × n_z
│ (1,1)   │ (1,2)   │ (1,3)   │  (n_z = number of altitudes)
├─────────┼─────────┼─────────┤
│ Block   │ Block   │ Block   │  Off-diagonal blocks (i1≠i2):
│ (2,1)   │ (2,2)   │ (2,3)   │  represent angular scattering (B matrix)
├─────────┼─────────┼─────────┤
│ Block   │ Block   │ Block   │  Diagonal blocks (i1=i2):
│ (3,1)   │ (3,2)   │ (3,3)   │  transport + loss + diffusion
└─────────┴─────────┴─────────┘
```

In the context of solving `f(Ie) = Mlhs*Ie - Q = 0`, the matrix `Mlhs` is the Jacobian:
```
Jacobian = ∂f/∂Ie = Mlhs
```

# Arguments
- `Ie`: pre-allocated output array [m⁻² s⁻¹] (n_z * n_angle) to store results
- `h_atm`: altitude grid [km]
- `μ`: cosine of pitch angle grid
- `matrices::TransportMatrices`: container with
    - `A`: electron loss rate [s⁻¹]
    - `B`: scattering matrix [s⁻¹] (n_z × n_angle × n_angle)
    - `D`: pitch-angle diffusion coefficient [s⁻¹] (n_angle,)
    - `Q`: source term [m⁻² s⁻¹] at each time step
    - `Ddiffusion`: spatial diffusion matrix (n_z × n_z)
- `iE`: current energy index
- `Ie_top`: boundary condition at top [m⁻² s⁻¹]
- `cache`: Cache object storing Mlhs, mapping, KLU, and differentiation matrices
- `first_iteration`: whether this is the first call (creates structure) or subsequent (reuses structure)
"""
function steady_state_scheme_optimized!(Ie, h_atm, μ, matrices, iE, Ie_top, cache; first_iteration = false)
    n_z = length(h_atm)
    n_angle = length(μ)

    # Extract matrices from container
    A = matrices.A
    B = matrices.B
    D = @view(matrices.D[iE, :])  # Extract D slice for current energy
    Q_slice = @view(matrices.Q[:, :, iE])  # Extract Q slice for current energy
    Ddiffusion = matrices.Ddiffusion

    # Compute or retrieve differentiation matrices
    if first_iteration
        # Spatial differentiation matrices
        h4diffu = [h_atm[1] - (h_atm[2] - h_atm[1]) ; h_atm]
        h4diffd = [h_atm ; h_atm[end] + (h_atm[end] - h_atm[end-1])]
        Ddz_Up   = spdiagm(-1 => -1 ./ diff(h4diffu[2:end]),
                            0 =>  1 ./ diff(h4diffu[1:end]))
        Ddz_Down = spdiagm( 0 => -1 ./ diff(h4diffd[1:end]),
                            1 =>  1 ./ diff(h4diffd[1:end-1]))

        # First time: create sparsity pattern and mapping
        cache.Mlhs = create_steady_state_sparsity_pattern(n_z, n_angle, μ, D, Ddiffusion)
        cache.mapping = create_steady_state_nzval_mapping(cache.Mlhs, n_z, n_angle)
        cache.Ddz_Up = Ddz_Up
        cache.Ddz_Down = Ddz_Down
    else
        # Reuse stored differentiation matrices (they don't change if h_atm doesn't change)
        Ddz_Up = cache.Ddz_Up
        Ddz_Down = cache.Ddz_Down
    end

    # Update matrix values (fast, no allocations)
    update_steady_state_matrix!(cache.Mlhs, cache.mapping, A, B, D, Ddiffusion, Ddz_Up, Ddz_Down, μ, h_atm)

    # Update or create KLU factorization
    if first_iteration
        cache.KLU = klu(cache.Mlhs)
    else
        klu!(cache.KLU, cache.Mlhs)
    end

    # Set boundary conditions
    index_top_bottom = sort(vcat(1:length(h_atm):(length(μ)*length(h_atm)),
                            length(h_atm):length(h_atm):(length(μ)*length(h_atm))))

    I_top_bottom = (Ie_top * [0, 1]')'
    Q_local = copy(Q_slice)
    Q_local[index_top_bottom] = I_top_bottom[:]

    # Solve system
    Ie .= cache.KLU \ Q_local

    Ie[Ie .< 0] .= 0; # the fluxes should never be negative

    return nothing
end
