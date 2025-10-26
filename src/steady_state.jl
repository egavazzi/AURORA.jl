using KLU: klu, klu!
using LinearAlgebra: Diagonal
using SparseArrays: spdiagm, sparse!, dropzeros!, findnz

function steady_state_scheme(h_atm, μ, matrices, iE, Ie_top, cache; first_iteration = false)
    Ie = Array{Float64}(undef, length(h_atm) * length(μ))

    # Extract matrices from container
    A = matrices.A
    B = matrices.B
    D = @view(matrices.D[iE, :])  # Extract D slice for current energy
    Q_slice = @view(matrices.Q[:, 1, iE])  # Extract Q slice for current energy (time index 1)
    Ddiffusion = matrices.Ddiffusion

    # Spatial differentiation matrices, for up and down streams
    # Here we tuck on a fictious height one tick below lowest height and on tick above highest
    # height just to make it possible to use the diff function.
    h4diffu = [h_atm[1] - (h_atm[2] - h_atm[1]) ; h_atm]
    h4diffd = [h_atm ; h_atm[end] + (h_atm[end] - h_atm[end-1])]
    Ddz_Up   = spdiagm(-1 => -1 ./ diff(h4diffu[2:end]),
                        0 =>  1 ./ diff(h4diffu[1:end]))
    Ddz_Down = spdiagm( 0 => -1 ./ diff(h4diffd[1:end]),
                        1 =>  1 ./ diff(h4diffd[1:end-1]))



    # Building the CN matrices
    Nz = length(h_atm)
    row_l = Vector{Int64}() # maybe using sizehint could help?
    col_l = Vector{Int64}()
    val_l = Vector{Float64}()
    for i1 in axes(B, 2)
        for i2 in axes(B, 2)
            A_tmp = A
            B_tmp = B[:, i1, i2]
            if i1 != i2
                tmp_lhs = -B_tmp
                tmp_lhs[1] = tmp_lhs[end] = 0

                idx_row = (i1 - 1) * Nz .+ (1:Nz)
                idx_col = (i2 - 1) * Nz .+ (1:Nz)
                append!(row_l, idx_row)
                append!(col_l, idx_col)
                append!(val_l, tmp_lhs)
            else
                if μ[i1] < 0    # downward fluxes
                    tmp_lhs =   μ[i1] .* Ddz_Down .+ Diagonal(A_tmp) .- D[i1] .* Ddiffusion .+ Diagonal(-B_tmp)
                    tmp_lhs[[1, end], :] .= 0
                    tmp_lhs[end, end] = 1
                else            # upward fluxes
                    tmp_lhs =   μ[i1] .* Ddz_Up .+ Diagonal(A_tmp) .- D[i1] .* Ddiffusion .+ Diagonal(-B_tmp)
                    tmp_lhs[[1, end], :] .= 0
                    tmp_lhs[end, end-1:end] = [-1, 1]
                end
                tmp_lhs[1, 1] = 1

                idx_row = (i1 - 1) * Nz .+ findnz(tmp_lhs)[1]
                idx_col = (i2 - 1) * Nz .+ findnz(tmp_lhs)[2]
                append!(row_l, idx_row)
                append!(col_l, idx_col)
                append!(val_l, findnz(tmp_lhs)[3])
            end
        end
    end
    Mlhs = sparse!(row_l, col_l, val_l)
    dropzeros!(Mlhs)    # for performance

    index_top_bottom = sort(vcat(1:length(h_atm):(length(μ)*length(h_atm)),
                            length(h_atm):length(h_atm):(length(μ)*length(h_atm))))

    if first_iteration
        cache.KLU = klu(Mlhs)
    else
        klu!(cache.KLU, Mlhs)
    end

    I_top_bottom = (Ie_top * [0, 1]')'
    Q_local = copy(Q_slice)
    Q_local[index_top_bottom] = I_top_bottom[:]
    Ie = cache.KLU \ Q_local

    Ie[Ie .< 0] .= 0; # the fluxes should never be negative

    return Ie
end
