using KLU: klu, klu!
using LinearAlgebra: Diagonal, ldiv!, mul!
using SparseArrays: spdiagm, sparse, dropzeros!, findnz

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


function Crank_Nicolson(t, h_atm, μ, v, matrices, iE, Ie_top, I0, cache; first_iteration = false)
    Ie = Array{Float64}(undef, length(h_atm) * length(μ), length(t))

    # Extract matrices from container
    A = matrices.A
    B = matrices.B
    D = @view(matrices.D[iE, :])  # Extract D slice for current energy
    Q_slice = @view(matrices.Q[:, :, iE])  # Extract Q slice for current energy
    Ddiffusion = matrices.Ddiffusion

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
    dt = t[2] - t[1]
    Ddt = Diagonal([1 ./ (v * dt) for i in h_atm])

    # Building the CN matrices
    Nz = length(h_atm)
    row_l = Vector{Int64}() # maybe using sizehint could help?
    col_l = Vector{Int64}()
    val_l = Vector{Float64}()
    row_r = Vector{Int64}()
    col_r = Vector{Int64}()
    val_r = Vector{Float64}()
    for i1 in axes(B, 2)
        for i2 in axes(B, 2)
            A_tmp = A
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
                    tmp_lhs =   μ[i1] .* Ddz_Down .+ Ddt .+ Diagonal(A_tmp/2) .- D[i1] .* Ddiffusion .+ Diagonal(-B_tmp/2)
                    tmp_rhs = - μ[i1] .* Ddz_Down .+ Ddt .- Diagonal(A_tmp/2) .+ D[i1] .* Ddiffusion .+ Diagonal( B_tmp/2)
                    tmp_lhs[[1, end], :] .= 0
                    tmp_lhs[end, end] = 1
                else            # upward fluxes
                    tmp_lhs =   μ[i1] .* Ddz_Up .+ Ddt .+ Diagonal(A_tmp/2) .- D[i1] .* Ddiffusion .+ Diagonal(-B_tmp/2)
                    tmp_rhs = - μ[i1] .* Ddz_Up .+ Ddt .- Diagonal(A_tmp/2) .+ D[i1] .* Ddiffusion .+ Diagonal( B_tmp/2)
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
    Mlhs = sparse!(row_l, col_l, val_l)
    Mrhs = sparse!(row_r, col_r, val_r)
    dropzeros!(Mlhs)    # for performance
    dropzeros!(Mrhs)    # for performance

    index_top_bottom = sort(vcat(1:length(h_atm):(length(μ)*length(h_atm)),
                            length(h_atm):length(h_atm):(length(μ)*length(h_atm))))

    i_t = 1
    Ie[:, 1] = I0
    Ie_finer = I0
    b = similar(Ie_finer)

    if first_iteration
        cache.KLU = klu(Mlhs)
    else
        klu!(cache.KLU, Mlhs)
    end

    for i_t in 1:length(t) - 1
        I_top_bottom = (@view(Ie_top[:, i_t]) * [0, 1]')'
        Q_local = (@view(Q_slice[:, i_t]) .+ @view(Q_slice[:, i_t + 1])) ./ 2
        Q_local[index_top_bottom] = I_top_bottom[:]

        mul!(b, Mrhs, Ie_finer)
        Ie_finer .= b
        Ie_finer .+= Q_local
        ldiv!(cache.KLU, Ie_finer)

        Ie[:, i_t + 1] = Ie_finer
    end
    Ie[Ie .< 0] .= 0; # the fluxes should never be negative

    return Ie
end
