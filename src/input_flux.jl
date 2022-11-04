using MAT

function Ie_top_from_file(t, E, n_loop, μ_center, filename)
    Ie_top = Array{Float64}(undef, length(μ_center), (n_loop - 1) * (length(t) - 1) + length(t), length(E))

    file = matopen(filename)
    Ie_top_raw = read(file, "Ie_total")
    close(file)

    # for constant input flux (e.g. first run), we need to resize the matrix from
    # [1, n_μ * [n_E, 1]] to [1, n_μ * [n_E, n_t]]
    if size(Ie_top[1], 2) == 1  
        for i_μ in eachindex(μ_center)
            Ie_top_raw[i_μ] = Ie_top_raw[i_μ] * ones(1, length(t) + (n_loop - 1) * (length(t) - 1)) # (e-/m²/s)
        end
    end

    # then we resize the matrix from [1, n_μ * [n_E, n_t]] to [n_μ, n_t, n_E] to be consistent with
    # the other flux matrices
    for i_μ in eachindex(μ_center)
        Ie_top[i_μ, :,  :] = Ie_top_raw[i_μ][1:length(E), :]' # (e-/m²/s)
    end

    return Ie_top
end


function Ie_top_flickering(t, E, dE, n_loop, μ_center, h_atm, BeamWeight_discrete, IeE_tot, z₀, E_min, f, Beams, modulation)
    Ie_top = zeros(length(μ_center), (n_loop - 1) * (length(t) - 1) + length(t), length(E))
    qₑ = 1.602176620898e-19
    i_Emin = findmin(abs.(E .- E_min))[2]   # find the index for the lower limit of the FAB
    
    z = z₀ * 1e3 - h_atm[end]   # distance between the source and the top of the ionosphere
    t_tot = t[1]:Float64(t.step):(t[end] * n_loop)
    t_shift₀ = z ./ (abs.(μ_center[Beams[1]]) * v_of_E(E[end]))

    for i_μ in eachindex(μ_center[Beams])
        for iE in eachindex(E)
            t_shift = z ./ (abs.(μ_center[Beams[i_μ]]) .* v_of_E.(E[iE]))
            t_shifted = t_tot .- (t_shift .- t_shift₀)
    
            IePnL = IeE_tot ./ qₑ ./ 
                    sum((E[i_Emin:end] .+ dE[i_Emin:end] ./ 2) .* dE[i_Emin:end]) .* dE
    
                if modulation == "sinus"
                    Ie_top[i_μ, :, iE] = IePnL[iE] .* (1 .- cos.(2*π*f/2 .* t_shifted) .^ 2) .* 
                                                        (t_shifted .> 0) .* (E[iE] > E[i_Emin]) .*
                                                        BeamWeight_discrete[i_μ] ./ sum(BeamWeight_discrete[Beams])
                elseif modulation == "square"
                    Ie_top[i_μ, :, iE] = IePnL[iE] .* (1 .+ square.(2*π*f .* t_shifted .- π/2) ./ 2) .* 
                                                        (t_shifted .> 0) .* (E[iE] > E[i_Emin]) .*
                                                        BeamWeight_discrete[i_μ] ./ sum(BeamWeight_discrete[Beams])
                end
        end
    end

    return Ie_top
end