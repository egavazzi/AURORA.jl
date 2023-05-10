using MAT

function Ie_top_from_old_matlab_file(t, E, n_loop, μ_center, filename)
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


function Ie_top_from_file(t, E, μ_center, n_loop, filename)
    Nt = (n_loop - 1) * (length(t) - 1) + length(t)

    ## load the file
    file = matopen(filename)
        Ie_top_raw = read(file, "Ie_total")
        t_top = read(file, "t")
    close(file)

    ## check that Ie_top is matching our simulation grid
    if size(Ie_top_raw, 1) != length(μ_center)
        error("""The incoming flux Ie_top is wrongly dimensioned. Check the θ dimension.\
        Remember that Ie_top should have the shape [n_μ, n_t, n_E].""")
    end
    if size(Ie_top_raw, 3) < length(E)
        error("""The incoming flux Ie_top is wrongly dimensioned. Check the E dimension.\
        Remember that Ie_top should have the shape [n_μ, n_t, n_E].""")
        error_grid = 1
    end

    ## check the time grid of Ie_top
    if length(t_top) == 1
        # we assume constant input flux and resize the matrix from
        # [n_μ, 1, n_E] to [n_μ, n_t, n_E]
        Ie_top = repeat(Ie_top_raw, outer=(1, Nt, 1))[:, :, 1:length(E)]
    else
        dt_top = t_top[2] - t_top[1]
        dt = t[2] - t[1]
        if dt_top > dt
            # the resolution in Ie_top is coarser than our simulation, we need to repeat elements
            dt_factor = dt_top / dt
            if !(dt_factor ≈ round(dt_factor))
                error("""Problem with the time resolution. The ratio of the dt of the \
                incoming flux over the dt of the simulation is not an integer. \n
                dt_incoming = $dt_top \n
                dt_simulation = $dt \n
                ratio = $dt_factor""")
            end
            dt_factor = round(Int, dt_factor)
            Ie_top_raw = repeat(Ie_top_raw, inner=(1, dt_factor, 1))
        elseif dt_top < dt
            # the resolution in Ie_top is finer than our simulation, we need to drop elements
            dt_factor = dt / dt_top
            if !(dt_factor ≈ round(dt_factor))
                error("""Problem with the time resolution. The ratio of the dt of the \
                simulation over the dt of the incoming flux is not an integer. \n
                dt_incoming = $dt_top \n
                dt_simulation = $dt \n
                ratio = $dt_factor""")
            end
            dt_factor = round(Int, dt_factor)
            Ie_top_raw = Ie_top_raw[:, 1:dt_factor:end, :]
        end

        if size(Ie_top_raw, 2) > Nt
            # in that case the array is too long and we need to cut it
            Ie_top_raw = Ie_top_raw[:, 1:Nt, :]
        elseif size(Ie_top_raw, 2) < Nt
            # in that case the array is too short and we fill it with zeros
            missing_data = zeros(length(μ_center), Nt - size(Ie_top_raw, 2), size(Ie_top_raw, 3))
            Ie_top_raw = cat(Ie_top_raw, missing_data; dims=2)
        end
        Ie_top = Ie_top_raw[:, :, 1:length(E)]
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



# Not sure about that one...
function Ie_top_Gaussian(t, E, n_loop, μ_center, h_atm, BeamWeight_discrete, I0, z₀, E0, dE0, Beams)
    Ie_top = zeros(length(μ_center), (n_loop - 1) * (length(t) - 1) + length(t), length(E))

    z = z₀ * 1e3 - h_atm[end]   # distance between the source and the top of the ionosphere
    t_tot = t[1]:Float64(t.step):(t[end] * n_loop)
    t_shift₀ = z ./ (abs.(μ_center[Beams[1]]) * v_of_E(E[end]))

    for i_μ in eachindex(μ_center[Beams])
        for iE in eachindex(E)
            t_shift = z ./ (abs.(μ_center[Beams[i_μ]]) .* v_of_E.(E[iE]))
            t_shifted = t_tot .- (t_shift .- t_shift₀)

            Ie_top[i_μ, :, iE] = I0 * exp(- (E[iE] - E0)^2 / dE0^2) .*
                                    (t_shifted .> 0) .*
                                    BeamWeight_discrete[i_μ] ./ sum(BeamWeight_discrete[Beams])
        end
    end

    return Ie_top
end



function Ie_top_constant(t, E, dE, n_loop, μ_center, h_atm, BeamWeight_discrete, IeE_tot, z₀, E_min, Beams, t0, t1)
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

            Ie_top[i_μ, :, iE] = IePnL[iE] .* f_smooth_transition.(t_shifted, t0, t1) .*
                                                (E[iE] > E[i_Emin]) .*
                                                BeamWeight_discrete[i_μ] ./ sum(BeamWeight_discrete[Beams])

        end
    end

    return Ie_top
end

function Ie_top_constant2(t, E, dE, n_loop, μ_center, h_atm, BeamWeight_discrete, IeE_tot, z₀, E_min, Beams, t0, t1)
    Ie_top = [zeros(length(μ_center), (n_loop - 1) * (length(t[i]) - 1) + length(t[i])) for i in eachindex(E)]
    qₑ = 1.602176620898e-19
    i_Emin = findmin(abs.(E .- E_min))[2]   # find the index for the lower limit of the FAB

    z = z₀ * 1e3 - h_atm[end]   # distance between the source and the top of the ionosphere
    t_shift₀ = z ./ (abs.(μ_center[Beams[1]]) * v_of_E(E[end]))

    for i_μ in eachindex(μ_center[Beams])
        for iE in eachindex(E)
            t_tot = t[iE][1]:Float64(t[iE].step):(t[iE][end] * n_loop)
            t_shift = z ./ (abs.(μ_center[Beams[i_μ]]) .* v_of_E.(E[iE]))
            t_shifted = t_tot .- (t_shift .- t_shift₀)

            IePnL = IeE_tot ./ qₑ ./
                    sum((E[i_Emin:end] .+ dE[i_Emin:end] ./ 2) .* dE[i_Emin:end]) .* dE

            Ie_top[iE][i_μ, :] = IePnL[iE] .* f_smooth_transition.(t_shifted, t0, t1) .*
                                                (E[iE] > E[i_Emin]) .*
                                                BeamWeight_discrete[i_μ] ./ sum(BeamWeight_discrete[Beams])

        end
    end

    return Ie_top
end
