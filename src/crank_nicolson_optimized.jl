# Optimized Crank-Nicolson scheme using direct nzval modification
# This version avoids allocations by reusing sparse matrix structures for both Mlhs and Mrhs

using KLU: klu, klu!
using LinearAlgebra: Diagonal, ldiv!, mul!
using SparseArrays: spdiagm, sparse, dropzeros!

# Note: d2M function is defined in crank_nicolson.jl and used here

"""
    create_crank_nicolson_sparsity_patterns(n_z, n_angle, μ, D, Ddiffusion)

Create the sparsity patterns for both Mlhs and Mrhs matrices.
The structure is the same for both, only values differ.

Returns (Mlhs, Mrhs) with correct sparsity structure.
"""
function create_crank_nicolson_sparsity_patterns(n_z, n_angle, μ, D, Ddiffusion)
    row_l = Vector{Int}()
    col_l = Vector{Int}()
    val_l = Vector{Float64}()
    row_r = Vector{Int}()
    col_r = Vector{Int}()
    val_r = Vector{Float64}()

    # Pre-allocate with estimated size
    max_nnz = n_angle * n_angle * 3 * n_z
    sizehint!(row_l, max_nnz)
    sizehint!(col_l, max_nnz)
    sizehint!(val_l, max_nnz)
    sizehint!(row_r, max_nnz)
    sizehint!(col_r, max_nnz)
    sizehint!(val_r, max_nnz)

    for i1 in 1:n_angle
        for i2 in 1:n_angle
            if i1 != i2
                # Off-diagonal blocks: diagonal elements only
                offset_row = (i1 - 1) * n_z
                offset_col = (i2 - 1) * n_z
                for i in 2:(n_z - 1)
                    # Mlhs
                    push!(row_l, offset_row + i)
                    push!(col_l, offset_col + i)
                    push!(val_l, 0.0)
                    # Mrhs
                    push!(row_r, offset_row + i)
                    push!(col_r, offset_col + i)
                    push!(val_r, 0.0)
                end
            else
                offset = (i1 - 1) * n_z

                # First row boundary condition
                push!(row_l, offset + 1)
                push!(col_l, offset + 1)
                push!(val_l, 1.0)
                # Mrhs boundary rows are zero (no entries needed in sparse matrix)

                if μ[i1] < 0    # downward fluxes
                    # Main tridiagonal part (skip first and last)
                    for i in 2:(n_z - 1)
                        # Diagonal
                        push!(row_l, offset + i)
                        push!(col_l, offset + i)
                        push!(val_l, 0.0)
                        push!(row_r, offset + i)
                        push!(col_r, offset + i)
                        push!(val_r, 0.0)

                        # Super-diagonal
                        push!(row_l, offset + i)
                        push!(col_l, offset + i + 1)
                        push!(val_l, 0.0)
                        push!(row_r, offset + i)
                        push!(col_r, offset + i + 1)
                        push!(val_r, 0.0)

                        # Sub-diagonal (from diffusion if present)
                        if D[i1] != 0 && Ddiffusion[i, i - 1] != 0
                            push!(row_l, offset + i)
                            push!(col_l, offset + i - 1)
                            push!(val_l, 0.0)
                            push!(row_r, offset + i)
                            push!(col_r, offset + i - 1)
                            push!(val_r, 0.0)
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
                        push!(row_r, offset + i)
                        push!(col_r, offset + i)
                        push!(val_r, 0.0)

                        # Sub-diagonal
                        push!(row_l, offset + i)
                        push!(col_l, offset + i - 1)
                        push!(val_l, 0.0)
                        push!(row_r, offset + i)
                        push!(col_r, offset + i - 1)
                        push!(val_r, 0.0)

                        # Super-diagonal (from diffusion if present)
                        if D[i1] != 0 && Ddiffusion[i, i + 1] != 0
                            push!(row_l, offset + i)
                            push!(col_l, offset + i + 1)
                            push!(val_l, 0.0)
                            push!(row_r, offset + i)
                            push!(col_r, offset + i + 1)
                            push!(val_r, 0.0)
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

    Mlhs = sparse(row_l, col_l, val_l, n_z * n_angle, n_z * n_angle)
    Mrhs = sparse(row_r, col_r, val_r, n_z * n_angle, n_z * n_angle)

    return Mlhs, Mrhs
end

"""
    create_crank_nicolson_nzval_mappings(Mlhs, Mrhs, n_z, n_angle)

Create mappings for both Mlhs and Mrhs matrices.
Returns (mapping_lhs, mapping_rhs).
"""
function create_crank_nicolson_nzval_mappings(Mlhs, Mrhs, n_z, n_angle)
    # Create mapping for Mlhs
    mapping_lhs = Matrix{Dict{Tuple{Symbol,Int},Int}}(undef, n_angle, n_angle)

    # Build reverse lookup for Mlhs
    nz_lookup_lhs = Dict{Tuple{Int,Int}, Int}()
    for col in 1:size(Mlhs, 2)
        for idx in Mlhs.colptr[col]:(Mlhs.colptr[col+1]-1)
            row = Mlhs.rowval[idx]
            nz_lookup_lhs[(row, col)] = idx
        end
    end

    # Create structured mapping for Mlhs
    for i1 in 1:n_angle
        for i2 in 1:n_angle
            block_map = Dict{Tuple{Symbol,Int},Int}()
            offset_row = (i1 - 1) * n_z
            offset_col = (i2 - 1) * n_z

            if i1 != i2
                # Off-diagonal: only diagonal elements
                for i in 2:(n_z - 1)
                    idx = get(nz_lookup_lhs, (offset_row + i, offset_col + i), nothing)
                    if idx !== nothing
                        block_map[(:diag, i)] = idx
                    end
                end
            else
                # Diagonal blocks
                offset = (i1 - 1) * n_z

                # Boundary conditions
                idx = get(nz_lookup_lhs, (offset + 1, offset + 1), nothing)
                if idx !== nothing
                    block_map[(:bc_first, 1)] = idx
                end

                # Interior points
                for i in 2:(n_z - 1)
                    idx = get(nz_lookup_lhs, (offset + i, offset + i), nothing)
                    if idx !== nothing
                        block_map[(:diag, i)] = idx
                    end

                    idx = get(nz_lookup_lhs, (offset + i, offset + i + 1), nothing)
                    if idx !== nothing
                        block_map[(:super, i)] = idx
                    end

                    idx = get(nz_lookup_lhs, (offset + i, offset + i - 1), nothing)
                    if idx !== nothing
                        block_map[(:sub, i)] = idx
                    end
                end

                # Last row boundary
                idx = get(nz_lookup_lhs, (offset + n_z, offset + n_z), nothing)
                if idx !== nothing
                    block_map[(:bc_last, n_z)] = idx
                end
                idx = get(nz_lookup_lhs, (offset + n_z, offset + n_z - 1), nothing)
                if idx !== nothing
                    block_map[(:bc_last_sub, n_z)] = idx
                end
            end

            mapping_lhs[i1, i2] = block_map
        end
    end

    # Create mapping for Mrhs (similar structure)
    mapping_rhs = Matrix{Dict{Tuple{Symbol,Int},Int}}(undef, n_angle, n_angle)

    # Build reverse lookup for Mrhs
    nz_lookup_rhs = Dict{Tuple{Int,Int}, Int}()
    for col in 1:size(Mrhs, 2)
        for idx in Mrhs.colptr[col]:(Mrhs.colptr[col+1]-1)
            row = Mrhs.rowval[idx]
            nz_lookup_rhs[(row, col)] = idx
        end
    end

    # Create structured mapping for Mrhs
    for i1 in 1:n_angle
        for i2 in 1:n_angle
            block_map = Dict{Tuple{Symbol,Int},Int}()
            offset_row = (i1 - 1) * n_z
            offset_col = (i2 - 1) * n_z

            if i1 != i2
                # Off-diagonal: only diagonal elements
                for i in 2:(n_z - 1)
                    idx = get(nz_lookup_rhs, (offset_row + i, offset_col + i), nothing)
                    if idx !== nothing
                        block_map[(:diag, i)] = idx
                    end
                end
            else
                # Diagonal blocks (no boundary entries in RHS)
                offset = (i1 - 1) * n_z

                # Interior points only
                for i in 2:(n_z - 1)
                    idx = get(nz_lookup_rhs, (offset + i, offset + i), nothing)
                    if idx !== nothing
                        block_map[(:diag, i)] = idx
                    end

                    idx = get(nz_lookup_rhs, (offset + i, offset + i + 1), nothing)
                    if idx !== nothing
                        block_map[(:super, i)] = idx
                    end

                    idx = get(nz_lookup_rhs, (offset + i, offset + i - 1), nothing)
                    if idx !== nothing
                        block_map[(:sub, i)] = idx
                    end
                end
            end

            mapping_rhs[i1, i2] = block_map
        end
    end

    return mapping_lhs, mapping_rhs
end

"""
    update_crank_nicolson_matrices!(Mlhs, Mrhs, mapping_lhs, mapping_rhs,
                                    A, B, D, Ddt, Ddiffusion, Ddz_Up, Ddz_Down, μ, h_atm)

Update both Mlhs and Mrhs using pre-computed mappings.
This is the fast path with zero allocations.
"""
function update_crank_nicolson_matrices!(Mlhs, Mrhs, mapping_lhs, mapping_rhs,
                                         A, B, D, Ddt, Ddiffusion, Ddz_Up, Ddz_Down, μ, h_atm)
    n_z = length(h_atm)
    n_angle = length(μ)
    nzval_lhs = Mlhs.nzval
    nzval_rhs = Mrhs.nzval

    for i1 in 1:n_angle
        for i2 in 1:n_angle
            B_tmp = @view B[:, i1, i2]
            block_map_lhs = mapping_lhs[i1, i2]
            block_map_rhs = mapping_rhs[i1, i2]

            if i1 != i2
                # Off-diagonal blocks
                for i in 2:(n_z - 1)
                    # Mlhs: -B/2
                    idx = get(block_map_lhs, (:diag, i), nothing)
                    if idx !== nothing
                        nzval_lhs[idx] = -B_tmp[i] / 2
                    end
                    # Mrhs: +B/2
                    idx = get(block_map_rhs, (:diag, i), nothing)
                    if idx !== nothing
                        nzval_rhs[idx] = B_tmp[i] / 2
                    end
                end
            else
                # Boundary conditions
                idx = get(block_map_lhs, (:bc_first, 1), nothing)
                if idx !== nothing
                    nzval_lhs[idx] = 1.0
                end

                if μ[i1] < 0    # downward fluxes
                    for i in 2:(n_z - 1)
                        # Mlhs diagonal: μ*Ddz_Down + Ddt + A/2 - D*Ddiffusion - B/2
                        val = μ[i1] * Ddz_Down[i, i] + Ddt[i, i] + A[i]/2 - B_tmp[i]/2
                        if D[i1] != 0 && Ddiffusion[i, i] != 0
                            val -= D[i1] * Ddiffusion[i, i]
                        end
                        idx = get(block_map_lhs, (:diag, i), nothing)
                        if idx !== nothing
                            nzval_lhs[idx] = val
                        end

                        # Mrhs diagonal: -μ*Ddz_Down + Ddt - A/2 + D*Ddiffusion + B/2
                        val = -μ[i1] * Ddz_Down[i, i] + Ddt[i, i] - A[i]/2 + B_tmp[i]/2
                        if D[i1] != 0 && Ddiffusion[i, i] != 0
                            val += D[i1] * Ddiffusion[i, i]
                        end
                        idx = get(block_map_rhs, (:diag, i), nothing)
                        if idx !== nothing
                            nzval_rhs[idx] = val
                        end

                        # Super-diagonal
                        val_lhs = μ[i1] * Ddz_Down[i, i + 1]
                        val_rhs = -μ[i1] * Ddz_Down[i, i + 1]
                        if D[i1] != 0 && Ddiffusion[i, i + 1] != 0
                            val_lhs -= D[i1] * Ddiffusion[i, i + 1]
                            val_rhs += D[i1] * Ddiffusion[i, i + 1]
                        end
                        idx = get(block_map_lhs, (:super, i), nothing)
                        if idx !== nothing
                            nzval_lhs[idx] = val_lhs
                        end
                        idx = get(block_map_rhs, (:super, i), nothing)
                        if idx !== nothing
                            nzval_rhs[idx] = val_rhs
                        end

                        # Sub-diagonal
                        if D[i1] != 0 && Ddiffusion[i, i - 1] != 0
                            val_lhs = -D[i1] * Ddiffusion[i, i - 1]
                            val_rhs = D[i1] * Ddiffusion[i, i - 1]
                            idx = get(block_map_lhs, (:sub, i), nothing)
                            if idx !== nothing
                                nzval_lhs[idx] = val_lhs
                            end
                            idx = get(block_map_rhs, (:sub, i), nothing)
                            if idx !== nothing
                                nzval_rhs[idx] = val_rhs
                            end
                        end
                    end

                    # Last row boundary
                    idx = get(block_map_lhs, (:bc_last, n_z), nothing)
                    if idx !== nothing
                        nzval_lhs[idx] = 1.0
                    end

                else  # upward fluxes
                    for i in 2:(n_z - 1)
                        # Mlhs diagonal: μ*Ddz_Up + Ddt + A/2 - D*Ddiffusion - B/2
                        val = μ[i1] * Ddz_Up[i, i] + Ddt[i, i] + A[i]/2 - B_tmp[i]/2
                        if D[i1] != 0 && Ddiffusion[i, i] != 0
                            val -= D[i1] * Ddiffusion[i, i]
                        end
                        idx = get(block_map_lhs, (:diag, i), nothing)
                        if idx !== nothing
                            nzval_lhs[idx] = val
                        end

                        # Mrhs diagonal: -μ*Ddz_Up + Ddt - A/2 + D*Ddiffusion + B/2
                        val = -μ[i1] * Ddz_Up[i, i] + Ddt[i, i] - A[i]/2 + B_tmp[i]/2
                        if D[i1] != 0 && Ddiffusion[i, i] != 0
                            val += D[i1] * Ddiffusion[i, i]
                        end
                        idx = get(block_map_rhs, (:diag, i), nothing)
                        if idx !== nothing
                            nzval_rhs[idx] = val
                        end

                        # Sub-diagonal
                        val_lhs = μ[i1] * Ddz_Up[i, i - 1]
                        val_rhs = -μ[i1] * Ddz_Up[i, i - 1]
                        if D[i1] != 0 && Ddiffusion[i, i - 1] != 0
                            val_lhs -= D[i1] * Ddiffusion[i, i - 1]
                            val_rhs += D[i1] * Ddiffusion[i, i - 1]
                        end
                        idx = get(block_map_lhs, (:sub, i), nothing)
                        if idx !== nothing
                            nzval_lhs[idx] = val_lhs
                        end
                        idx = get(block_map_rhs, (:sub, i), nothing)
                        if idx !== nothing
                            nzval_rhs[idx] = val_rhs
                        end

                        # Super-diagonal
                        if D[i1] != 0 && Ddiffusion[i, i + 1] != 0
                            val_lhs = -D[i1] * Ddiffusion[i, i + 1]
                            val_rhs = D[i1] * Ddiffusion[i, i + 1]
                            idx = get(block_map_lhs, (:super, i), nothing)
                            if idx !== nothing
                                nzval_lhs[idx] = val_lhs
                            end
                            idx = get(block_map_rhs, (:super, i), nothing)
                            if idx !== nothing
                                nzval_rhs[idx] = val_rhs
                            end
                        end
                    end

                    # Last row boundary
                    idx = get(block_map_lhs, (:bc_last_sub, n_z), nothing)
                    if idx !== nothing
                        nzval_lhs[idx] = -1.0
                    end
                    idx = get(block_map_lhs, (:bc_last, n_z), nothing)
                    if idx !== nothing
                        nzval_lhs[idx] = 1.0
                    end
                end
            end
        end
    end

    return Mlhs, Mrhs
end

"""
    Crank_Nicolson_optimized!(Ie, t, h_atm, μ, v, matrices, iE, Ie_top, I0, cache; first_iteration = false)

Optimized Crank-Nicolson time-stepping scheme using direct nzval modification.
This is an in-place version that modifies `Ie` directly to avoid allocations.

On first iteration, creates the sparsity pattern and mapping which are stored in cache.
On subsequent iterations, only updates the nzval array directly.

# Mathematical Background

The time-dependent electron transport equation is:
```
∂Ie/∂t + μ*v ∂Ie/∂z + A*Ie - ∫B*Ie'dΩ' - D*∇²Ie = Q
```

The Crank-Nicolson scheme uses implicit time-stepping with second-order accuracy:
```
(Ie^(n+1) - Ie^n)/Δt = [RHS^(n+1) + RHS^n] / 2
```

This leads to two matrices:
```
Mlhs * Ie^(n+1) = Mrhs * Ie^n + (Q^(n+1) + Q^n)/2
     ↑                  ↑
Implicit part      Explicit part
```

Where:
```
Mlhs = Ddt + μ*Ddz/2 + A/2 - B/2 - D*Ddiffusion/2
Mrhs = Ddt - μ*Ddz/2 - A/2 + B/2 + D*Ddiffusion/2
```

Both matrices have the same block structure as in steady-state:
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
- `Ie`: pre-allocated output array [m⁻² s⁻¹], size (n_z * n_angle × n_t) to store results
- `t`: time grid [s]
- `h_atm`: altitude grid [km]
- `μ`: cosine of pitch angle grid
- `v`: electron velocity [km/s]
- `matrices::TransportMatrices`: container with
    - `A`: electron loss rate [s⁻¹]
    - `B`: scattering matrix [s⁻¹], size (n_z × n_angle × n_angle)
    - `D`: pitch-angle diffusion coefficient [s⁻¹], size (n_angle,)
    - `Q`: source term [m⁻² s⁻¹] at each time step
    - `Ddiffusion`: spatial diffusion matrix (n_z × n_z)
- `iE`: current energy index
- `Ie_top`: boundary condition at top [m⁻² s⁻¹] at each time step
- `I0`: initial condition [m⁻² s⁻¹]
- `cache`: Cache object (must have fields for Mlhs, Mrhs, mappings, KLU, diff matrices)
- `first_iteration`: whether this is the first call

"""
function Crank_Nicolson_optimized!(Ie, t, h_atm, μ, v, matrices, iE, Ie_top, I0, cache; first_iteration = false)
    n_z = length(h_atm)
    n_angle = length(μ)

    # Extract matrices from container
    A = matrices.A
    B = matrices.B
    D = @view(matrices.D[iE, :])  # Extract D slice for current energy
    Q_slice = @view(matrices.Q[:, :, iE])  # Extract Q slice for current energy
    Ddiffusion = matrices.Ddiffusion

    # Temporal differentiation matrix
    dt = t[2] - t[1]
    Ddt = Diagonal([1 ./ (v * dt) for i in h_atm])

    if first_iteration
        # Spatial differentiation matrices
        h4diffu = [h_atm[1] - (h_atm[2] - h_atm[1]); h_atm]
        h4diffd = [h_atm; h_atm[end] + (h_atm[end] - h_atm[end - 1])]
        Ddz_Up = spdiagm(-1 => -1 ./ (2 .* diff(h4diffu[2:end])),
                         0 => 1 ./ (2 .* diff(h4diffu[1:end])))
        Ddz_Down = spdiagm(0 => -1 ./ (2 .* diff(h4diffd[1:end])),
                           1 => 1 ./ (2 .* diff(h4diffd[1:(end - 1)])))

        # First time: create sparsity patterns and mappings
        cache.Mlhs, cache.Mrhs = create_crank_nicolson_sparsity_patterns(n_z, n_angle, μ, D, Ddiffusion)
        cache.mapping_lhs, cache.mapping_rhs = create_crank_nicolson_nzval_mappings(cache.Mlhs, cache.Mrhs, n_z, n_angle)
        cache.Ddz_Up = Ddz_Up
        cache.Ddz_Down = Ddz_Down
    else
        # Reuse cached matrices
        Ddz_Up = cache.Ddz_Up
        Ddz_Down = cache.Ddz_Down
    end

    # Update matrix values (fast, no allocations)
    update_crank_nicolson_matrices!(cache.Mlhs, cache.Mrhs, cache.mapping_lhs, cache.mapping_rhs,
                                     A, B, D, Ddt, Ddiffusion, Ddz_Up, Ddz_Down, μ, h_atm)

    # Boundary indices
    index_top_bottom = sort(vcat(1:n_z:(n_angle*n_z),
                            n_z:n_z:(n_angle*n_z)))

    # Initialize
    Ie[:, 1] = I0
    Ie_finer = copy(I0)
    b = similar(Ie_finer)

    # Create or update KLU factorization
    if first_iteration
        cache.KLU = klu(cache.Mlhs)
    else
        klu!(cache.KLU, cache.Mlhs)
    end

    # Time-stepping loop
    for i_t in 1:(length(t) - 1)
        I_top_bottom = (@view(Ie_top[:, i_t]) * [0, 1]')'
        Q_local = (@view(Q_slice[:, i_t]) .+ @view(Q_slice[:, i_t + 1])) ./ 2
        Q_local[index_top_bottom] = I_top_bottom[:]

        # Crank-Nicolson step: Mlhs * Ie^(n+1) = Mrhs * Ie^n + Q
        mul!(b, cache.Mrhs, Ie_finer)
        Ie_finer .= b
        Ie_finer .+= Q_local
        ldiv!(cache.KLU, Ie_finer)

        Ie[:, i_t + 1] = Ie_finer
    end

    Ie[Ie .< 0] .= 0  # Fluxes should never be negative

    return nothing
end
