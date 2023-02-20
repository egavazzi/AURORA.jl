using LinearAlgebra
using SparseArrays

function d2M(z)
    dzd = z[2:end-1] - z[1:end-2]
    dzu = z[3:end]   - z[2:end-1]

    dsup  = [2 ./ (dzd .* (dzd + dzu)) ; 0]
    dMain = [0 ; -2 ./ (dzd .* dzu) ; 0]
    dsub  = [0 ; 2 ./ (dzu .* (dzd + dzu))]

    D2M = spdiagm( -1 => dsup,
                    0 => dMain,
                    1 => dsub)
    return D2M
end



using LinearSolve
function Crank_Nicolson(t, h_atm, μ, v, A, B, D, Q, Ie_top, I0)
    Ie = Array{Float64}(undef, length(h_atm) * length(μ), length(t))
    # Ie = zeros(length(h_atm) * length(μ), length(t))

    dt = t[2] - t[1]
    dz = h_atm[2] - h_atm[1] # smallest value of dz

    # The Courant-Freidrichs-Lewy (CFL) number normally hase to be small (<4) to ensure numerical
    # stability. However, as a Crank-Nicolson scheme is always stable, we can take a bigger CFL. We
    # should be careful about numerical accuracy though.
    # For Gaussian inputs (or similar), it seems that the CFL can be set to 64 without major effects
    # on the results, while reducing computational time tremendously
    CFL = v * dt / dz
    n_factors = 2 .^ collect(0:22)
    iFactor = 1
    # This while loop effectively reduces dt by a factor of 2 at each iteration and check if the new
    # CFL is < 64. If not, it continues reducing dt.
    t_finer = t
    while (CFL > 64) && (iFactor < length(n_factors))
        t_finer = range(t[1], t[end], length(t) * n_factors[iFactor] + 1 - n_factors[iFactor])
        dt = t_finer[2] - t_finer[1]
        CFL = v * dt / dz
        iFactor += 1
    end
    CFL_factor = n_factors[max(1, iFactor - 1)]

    # Spatial differentiation matrices, for up and down streams
    # Here we tuck on a fictious height one tick below lowest height and on tick above highest
    # height just to make it possible to use the diff function.
    h4diffu = [h_atm[1] - (h_atm[2] - h_atm[1]) ; h_atm]
    h4diffd = [h_atm ; h_atm[end] + (h_atm[end] - h_atm[end-1])]
    Ddz_Up   = spdiagm(-1 => -1 ./ (2 .* diff(h4diffu[2:end])),
                        0 =>  1 ./ (2 .* diff(h4diffu[1:end])))
    Ddz_Down = spdiagm( 0 => -1 ./ (2 .* diff(h4diffd[1:end])),
                        1 =>  1 ./ (2 .* diff(h4diffd[1:end-1])))

    # Temporal differentiation matrix
    Ddt = Diagonal([1 ./ (v * dt) for i in h_atm])

    # Diffusion operator
    Ddiffusion = d2M(h_atm)
    Ddiffusion[1, 1] = 0

    # Building the CN matrices
    Mlhs = Array{Float64}(undef, 0, length(μ)*length(h_atm))
    Mrhs = Array{Float64}(undef, 0, length(μ)*length(h_atm))
    for i1 in axes(B, 2)
        MBlhs = Array{Float64}(undef, length(h_atm), 0)
        MBrhs = Array{Float64}(undef, length(h_atm), 0)
        for i2 in axes(B, 2)
            B_tmp = B[:, i1, i2]
            if i1 != i2
                tmp_lhs = Diagonal(-B_tmp/2)
                tmp_lhs[[1, end], :] .= 0
                tmp_rhs = Diagonal( B_tmp/2)
                tmp_rhs[[1, end], :] .= 0

                MBlhs = sparse_hcat(MBlhs, tmp_lhs)
                MBrhs = sparse_hcat(MBrhs, tmp_rhs)
            else
                if μ[i1] < 0    # downward fluxes
                    tmp_lhs =   μ[i1] .* Ddz_Down .+ Ddt .+ Diagonal(A/2) .- D[i1] .* Ddiffusion .+ Diagonal(-B_tmp/2)
                    tmp_rhs = - μ[i1] .* Ddz_Down .+ Ddt .- Diagonal(A/2) .+ D[i1] .* Ddiffusion .+ Diagonal( B_tmp/2)
                    tmp_lhs[[1, end], :] .= 0
                    tmp_lhs[end, end] = 1
                else            # upward fluxes
                    tmp_lhs =   μ[i1] .* Ddz_Up .+ Ddt .+ Diagonal(A/2) .- D[i1] .* Ddiffusion .+ Diagonal(-B_tmp/2)
                    tmp_rhs = - μ[i1] .* Ddz_Up .+ Ddt .- Diagonal(A/2) .+ D[i1] .* Ddiffusion .+ Diagonal( B_tmp/2)
                    tmp_lhs[[1, end], :] .= 0
                    tmp_lhs[end, end-1:end] = [-1, 1]
                end
                tmp_rhs[[1, end], :] .= 0
                tmp_lhs[1, 1] = 1

                MBlhs = sparse_hcat(MBlhs, tmp_lhs)
                MBrhs = sparse_hcat(MBrhs, tmp_rhs)
            end
        end
        Mlhs = sparse_vcat(Mlhs, MBlhs)
        Mrhs = sparse_vcat(Mrhs, MBrhs)
    end

    index_top_bottom = sort(vcat(1:length(h_atm):(length(μ)*length(h_atm)),
                            length(h_atm):length(h_atm):(length(μ)*length(h_atm))))

    i_t = 1
    Ie[:, 1] = I0
    Ie_finer = I0
    b = similar(Ie_finer)
    prob = LinearProblem(Mlhs, Ie_finer)
    linsolve = init(prob)
    sol1 = solve(linsolve)



    for i_t_finer in 2:length(t_finer)
        I_top_bottom = (@view(Ie_top[:, i_t]) * [0, 1]')'
        Q_local = (@view(Q[:, i_t]) .+ @view(Q[:, i_t + 1])) ./ 2
        Q_local[index_top_bottom] = I_top_bottom[:]

        mul!(b, Mrhs, Ie_finer)
        b .+= Q_local

        linsolve = LinearSolve.set_b(sol1.cache, b)
        sol2 = solve(linsolve)

        Ie_finer = sol2.u

        if rem(i_t_finer, CFL_factor) == 1 || CFL_factor == 1
            i_t = i_t + 1
            Ie[:, i_t] = Ie_finer
        end
    end
    Ie[Ie .< 0] .= 0; # the fluxes should never be negative
    return Ie
end


using KLU
function Crank_Nicolson_Optimized(t, h_atm, μ, v, A, B, D, Q, Ie_top, I0)
    Ie = Array{Float64}(undef, length(h_atm) * length(μ), length(t))
    # Ie = zeros(length(h_atm) * length(μ), length(t))

    dt = t[2] - t[1]
    dz = h_atm[2] - h_atm[1] # smallest value of dz

    # The Courant-Freidrichs-Lewy (CFL) number normally hase to be small (<4) to ensure numerical
    # stability. However, as a Crank-Nicolson scheme is always stable, we can take a bigger CFL. We
    # should be careful about numerical accuracy though.
    # For Gaussian inputs (or similar), it seems that the CFL can be set to 64 without major effects
    # on the results, while reducing computational time tremendously
    CFL = v * dt / dz
    n_factors = 2 .^ collect(0:22)
    iFactor = 1
    # This while loop effectively reduces dt by a factor of 2 at each iteration and check if the new
    # CFL is < 64. If not, it continues reducing dt.
    t_finer = t
    while (CFL > 64) && (iFactor < length(n_factors))
        t_finer = range(t[1], t[end], length(t) * n_factors[iFactor] + 1 - n_factors[iFactor])
        dt = t_finer[2] - t_finer[1]
        CFL = v * dt / dz
        iFactor += 1
    end
    CFL_factor = n_factors[max(1, iFactor - 1)]

    # Spatial differentiation matrices, for up and down streams
    # Here we tuck on a fictious height one tick below lowest height and on tick above highest
    # height just to make it possible to use the diff function.
    h4diffu = [h_atm[1] - (h_atm[2] - h_atm[1]) ; h_atm]
    h4diffd = [h_atm ; h_atm[end] + (h_atm[end] - h_atm[end-1])]
    Ddz_Up   = spdiagm(-1 => -1 ./ (2 .* diff(h4diffu[2:end])),
                        0 =>  1 ./ (2 .* diff(h4diffu[1:end])))
    Ddz_Down = spdiagm( 0 => -1 ./ (2 .* diff(h4diffd[1:end])),
                        1 =>  1 ./ (2 .* diff(h4diffd[1:end-1])))

    # Temporal differentiation matrix
    Ddt = Diagonal([1 ./ (v * dt) for i in h_atm])

    # Diffusion operator
    Ddiffusion = d2M(h_atm)
    Ddiffusion[1, 1] = 0

    # Building the CN matrices
    Nz = length(h_atm)
    row_l = Vector{Int64}()
    col_l = Vector{Int64}()
    val_l = Vector{Float64}()
    row_r = Vector{Int64}()
    col_r = Vector{Int64}()
    val_r = Vector{Float64}()
    for i1 in axes(B, 2)
        for i2 in axes(B, 2)
            B_tmp = B[:, i1, i2]
            if i1 != i2
                tmp_lhs = -B_tmp/2
                tmp_lhs[1] = tmp_lhs[end] = 0

                tmp_rhs = B_tmp/2
                tmp_rhs[1] = tmp_rhs[end] = 0

                idx_row = (i1 - 1) * Nz .+ (1:Nz)
                idx_col = (i2 - 1) * Nz .+ (1:Nz)
                append!(row_l, idx_row)
                append!(col_l, idx_col)
                append!(val_l, tmp_lhs)
                append!(row_r, idx_row)
                append!(col_r, idx_col)
                append!(val_r, tmp_rhs)
            else
                if μ[i1] < 0    # downward fluxes
                    tmp_lhs =   μ[i1] .* Ddz_Down .+ Ddt .+ Diagonal(A/2) .- D[i1] .* Ddiffusion .+ Diagonal(-B_tmp/2)
                    tmp_rhs = - μ[i1] .* Ddz_Down .+ Ddt .- Diagonal(A/2) .+ D[i1] .* Ddiffusion .+ Diagonal( B_tmp/2)
                    tmp_lhs[[1, end], :] .= 0
                    tmp_lhs[end, end] = 1
                else            # upward fluxes
                    tmp_lhs =   μ[i1] .* Ddz_Up .+ Ddt .+ Diagonal(A/2) .- D[i1] .* Ddiffusion .+ Diagonal(-B_tmp/2)
                    tmp_rhs = - μ[i1] .* Ddz_Up .+ Ddt .- Diagonal(A/2) .+ D[i1] .* Ddiffusion .+ Diagonal( B_tmp/2)
                    tmp_lhs[[1, end], :] .= 0
                    tmp_lhs[end, end-1:end] = [-1, 1]
                end
                tmp_rhs[[1, end], :] .= 0
                tmp_lhs[1, 1] = 1

                idx_row = (i1 - 1) * Nz .+ findnz(tmp_lhs)[1]
                idx_col = (i2 - 1) * Nz .+ findnz(tmp_lhs)[2]
                append!(row_l, idx_row)
                append!(col_l, idx_col)
                append!(val_l, findnz(tmp_lhs)[3])
                idx_row = (i1 - 1) * Nz .+ findnz(tmp_rhs)[1]
                idx_col = (i2 - 1) * Nz .+ findnz(tmp_rhs)[2]
                append!(row_r, idx_row)
                append!(col_r, idx_col)
                append!(val_r, findnz(tmp_rhs)[3])
            end
        end
    end
    Mlhs = sparse(row_l, col_l, val_l)
    Mrhs = sparse(row_r, col_r, val_r)
    dropzeros!(Mlhs)    # for the performance of the next calculations
    dropzeros!(Mrhs)    # for the performance of the next calculations

    index_top_bottom = sort(vcat(1:length(h_atm):(length(μ)*length(h_atm)),
                            length(h_atm):length(h_atm):(length(μ)*length(h_atm))))

    i_t = 1
    Ie[:, 1] = I0
    Ie_finer = I0
    b = similar(Ie_finer)

    # klu!(AAA, Mlhs)
    AAA = klu(Mlhs)

    for i_t_finer in 2:length(t_finer)
        I_top_bottom = (@view(Ie_top[:, i_t]) * [0, 1]')'
        Q_local = (@view(Q[:, i_t]) .+ @view(Q[:, i_t + 1])) ./ 2
        Q_local[index_top_bottom] = I_top_bottom[:]

        mul!(b, Mrhs, Ie_finer)
        Ie_finer .= b
        Ie_finer .+= Q_local
        ldiv!(AAA, Ie_finer)

        if rem(i_t_finer, CFL_factor) == 1 || CFL_factor == 1
            i_t = i_t + 1
            Ie[:, i_t] = Ie_finer
        end
    end
    Ie[Ie .< 0] .= 0; # the fluxes should never be negative
    return Ie
end
