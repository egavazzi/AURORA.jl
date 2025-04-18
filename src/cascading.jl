using Dates: Dates, now
using HCubature: hcubature
using Interpolations: linear_interpolation
using MAT: matopen
using ProgressMeter: Progress, next!


# function is_cascading_loaded(Q, E4Q, E_secondary)
#     # check1: is Q empty? if yes we have to load
#     check1 = isempty(Q)

#     # check2: are E4Q and E_secondary the same? if not we have to load
#     check2 = E4Q[1:min(length(E_secondary), length(E4Q))] !=
#                 E_secondary[1:min(length(E_secondary), length(E4Q))]

#     # check3: is E_secondary larger than E4Q? if yes we have to load
#     # (note that the opposite (E4Q > E_secondary) is not a problem)
#     check3 = length(E_secondary) > length(E4Q)

#     return check1, check2, check3
# end


let Q = [], E4Q , E_ionizations
    global function cascading_N2(E_secondary, E_primary, E_ionization, s_or_c)

        # First we check if the cascading spectra has already been loaded by checking the
        # static/persistent local variables (namely Q, E4Q and E_ionizations)
        # check1, check2, check3 = is_cascading_loaded(Q, E4Q, E_secondary)

        E_hat = 11.4
        # if check1 || check2 || check3
        if isempty(Q) || length(E_secondary) > length(E4Q) || E4Q[1:length(E_secondary)] != E_secondary
            # then we try to find a cascading spectra file with matching energy grid E_secondary
            found_them = 0
            cascading_files = readdir(pkgdir(AURORA, "internal_data", "e_cascading", "N2"))
            for i1 in eachindex(cascading_files)
                if !isdir(cascading_files[i1])
                    try
                        filename = pkgdir(AURORA, "internal_data", "e_cascading", "N2", cascading_files[i1])
                        file = matopen(filename)
                        E4Q = read(file, "E4Q")
                        close(file)
                        if length(E_secondary) <= length(E4Q) && E4Q[1:length(E_secondary)] == E_secondary
                            println("Loading cascading-matrices from file: ", cascading_files[i1])
                            file = matopen(filename)
                            Q = read(file, "Q")
                            E_ionizations = read(file, "E_ionizations")
                            close(file)
                            found_them = 1
                        break
                        end
                    catch
                    end
                end
            end
            if found_them == 0
                println("Could not find a file with matching energy grid.")
                println("Starting to calculate the requested cascading-matrices.")

                # Precalculating the degrading spectra
                E_ionizations = [15.581, 16.73, 18.75, 24, 42]
                Q = zeros(length(E_secondary), length(E_secondary), length(E_ionizations))
                E4Q = E_secondary
                dE = diff(E_secondary); dE = dE[[1:end; end]]
                println("Pre-calculating all energy-degradations for e - N2-ionizing collisions.")
                println("Starting at ", Dates.format(now(), "HH:MM:SS"))
                for i1 in length(E_ionizations):-1:1
                    iLim = findfirst(x -> x > E_ionizations[i1], E_secondary)

                    # Lower boundary of primary electrons, This correctly accounts for the
                    # limit when Edeg[i3 + 1] is larger than Eprimary[i2 + 1] - Eionizations[i1],
                    # cutting a corner off the [Ei3 - (Ei3 + dEi3)] x [Ei2 + dEi2] square
                    # we're integrating over.
                    Emin(E) = max(Eprime[1], E + E_ionizations[i1])

                    # Change variable E[2] so that we can integrate 'fun' over an hypercube later
                    Eprime = Vector{Float64}    # we have to define an empty Eprime here so it can be
                                                # used in the function definition
                    var_hypercube(E1, E2) = Emin(E1) .+ E2 .* (Eprime[end] .- Emin(E1))
                    # Compute the jacobian associated to the change of variable above
                    jacobian(E1) = (Eprime[end] .- Emin(E1))

                    # Function fun... TODO
                    fun(E) = 1 ./
                            ((var_hypercube(E[1], E[2]) .- E_ionizations[i1] .- E[1]).^2 .+
                            E_hat^2) .* jacobian(E[1])

                    p = Progress(length(E_secondary) - iLim; desc=string("Doing level ", i1), color=:blue);
                    for i2 in iLim:length(E_secondary)
                        iHalf = findlast(x -> x < (E_secondary[i2] - E_ionizations[i1]) / 2, E_secondary)
                        if !isnothing(iHalf) # otherwise we go to next i2
                            Eprime = [E_secondary[i2], E_secondary[i2] + dE[i2]]

                            for i3 in iHalf:(i2 - 1)
                                # Correct upper Edeg integration limit, this cuts the region into a
                                # triangle when needed.
                                Edegmax = min(E_secondary[i3] + dE[i3], E_secondary[i2] + dE[i2] - E_ionizations[i1])
                                # This if-condition just makes sure that we have the physically meaningful
                                # limits, i.e. we integrate in increasing Edeg - which is not otherwise
                                # guaranteed.
                                if Edegmax > E_secondary[i3]
                                    Q[i2, i3, i1] = hcubature(fun, [E_secondary[i3], 0],
                                                                    [Edegmax, 1], rtol = 1e-4)[1]
                                end
                            end
                        end
                        next!(p)
                    end
                    # println("Done with level ", i1, " at: ", Dates.format(now(), "HH:MM:SS"))
                end

                # Save the results for future use
                filename = pkgdir(AURORA, "internal_data", "e_cascading", "N2",
                                    string("CascadingSpecN2ionization_",
                                    Dates.format(now(), "yyyymmdd-HHMMSS"),
                                    ".mat"))
                file = matopen(filename, "w")
                write(file, "Q", Q)
                write(file, "E4Q", E4Q)
                write(file, "E_ionizations", E_ionizations)
                close(file)
            end
        end


        if s_or_c == "s"
            # Calculating the spectra of the secondary e-
            secondary_e = 1 ./ (E_hat^2 .+ E_secondary.^2) .* (E_secondary .< (E_primary - E_ionization) / 2)

            return secondary_e
        elseif s_or_c == "c"
            # Calculating the spectra of the degrading primary e-
            iLevel = findmin(abs.(E_ionization .- E_ionizations))[2]
            iPrimary = findmin(abs.(E4Q .- E_primary))[2]
            primary_e = Q[iPrimary, 1:length(E_secondary), iLevel]

            return primary_e
        end
    end
end


let Q = [], E4Q , E_ionizations
    global function cascading_O2(E_secondary, E_primary, E_ionization, s_or_c)

        # First we check if the cascading spectra has already been loaded by checking the
        # static/persistent local variables (namely Q, E4Q and E_ionizations)
        # check1, check2, check3 = is_cascading_loaded(Q, E4Q, E_secondary)

        E_hat = 15.2
        # if check1 || check2 || check3
        if isempty(Q) || length(E_secondary) > length(E4Q) || E4Q[1:length(E_secondary)] != E_secondary
            # then we try to find a cascading spectra file with matching energy grid E_secondary
            found_them = 0
            cascading_files = readdir(pkgdir(AURORA, "internal_data", "e_cascading", "O2"))
            for i1 in eachindex(cascading_files)
                try
                    if !isdir(cascading_files[i1])
                        filename = pkgdir(AURORA, "internal_data", "e_cascading", "O2", cascading_files[i1])
                        file = matopen(filename)
                        E4Q = read(file, "E4Q")
                        close(file)
                        if length(E_secondary) <= length(E4Q) && E4Q[1:length(E_secondary)] == E_secondary
                            println("Loading cascading-matrices from file: ", cascading_files[i1])
                            file = matopen(filename)
                            Q = read(file, "Q")
                            E_ionizations = read(file, "E_ionizations")
                            close(file)
                            found_them = 1
                        break
                        end
                    end
                catch
                end
            end
            if found_them == 0
                println("Could not find a file with matching energy grid.")
                println("Starting to calculate the requested cascading-matrices.")

                # Precalculating the degrading spectra
                E_ionizations = [12.072, 16.1, 16.9, 18.2, 18.9, 32.51]
                Q = zeros(length(E_secondary), length(E_secondary), length(E_ionizations))
                E4Q = E_secondary
                dE = diff(E_secondary); dE = dE[[1:end; end]]
                println("Pre-calculating all energy-degradations for e - O2-ionizing collisions.")
                println("Starting at ", Dates.format(now(), "HH:MM:SS"))
                for i1 in length(E_ionizations):-1:1
                    iLim = findfirst(x -> x > E_ionizations[i1], E_secondary)

                    # Lower boundary of primary electrons, This correctly accounts for the
                    # limit when Edeg[i3 + 1] is larger than Eprimary[i2 + 1] - Eionizations[i1],
                    # cutting a corner off the [Ei3 - (Ei3 + dEi3)] x [Ei2 + dEi2] square
                    # we're integrating over.
                    Emin(E) = max(Eprime[1], E + E_ionizations[i1])

                    # Change variable E[2] so that we can integrate 'fun' over an hypercube later
                    Eprime = Vector{Float64}    # we have to define an empty Eprime here so it can be
                                                # used in the function definition
                    var_hypercube(E1, E2) = Emin(E1) .+ E2 .* (Eprime[end] .- Emin(E1))
                    # Compute the jacobian associated to the change of variable above
                    jacobian(E1) = (Eprime[end] .- Emin(E1))

                    # Function fun... TODO
                    fun(E) = 1 ./
                            ((var_hypercube(E[1], E[2]) .- E_ionizations[i1] .- E[1]).^2 .+
                            E_hat^2) .* jacobian(E[1])

                    p = Progress(length(E_secondary) - iLim; desc=string("Doing level ", i1), color=:blue);
                    for i2 in iLim:length(E_secondary)
                        iHalf = findlast(x -> x < (E_secondary[i2] - E_ionizations[i1]) / 2, E_secondary)
                        if !isnothing(iHalf) # otherwise we go to next i2
                            Eprime = [E_secondary[i2], E_secondary[i2] + dE[i2]]

                            for i3 in iHalf:(i2 - 1)
                                # Correct upper Edeg integration limit, this cuts the region into a
                                # triangle when needed.
                                Edegmax = min(E_secondary[i3] + dE[i3], E_secondary[i2] + dE[i2] - E_ionizations[i1])
                                # This if-condition just makes sure that we have the physically meaningful
                                # limits, i.e. we integrate in increasing Edeg - which is not otherwise
                                # guaranteed.
                                if Edegmax > E_secondary[i3]
                                    Q[i2, i3, i1] = hcubature(fun, [E_secondary[i3], 0],
                                                                    [Edegmax, 1], rtol = 1e-4)[1]
                                end
                            end
                        end
                        next!(p)
                    end
                    # println("Done with level ", i1, " at: ", Dates.format(now(), "HH:MM:SS"))
                end

                # Save the results for future use
                filename = pkgdir(AURORA, "internal_data", "e_cascading", "O2",
                                    string("CascadingSpecO2ionization_",
                                    Dates.format(now(), "yyyymmdd-HHMMSS"),
                                    ".mat"))
                file = matopen(filename, "w")
                write(file, "Q", Q)
                write(file, "E4Q", E4Q)
                write(file, "E_ionizations", E_ionizations)
                close(file)
            end
        end



        if s_or_c == "s"
            # Calculating the spectra of the secondary e-
            secondary_e = 1 ./ (E_hat^2 .+ E_secondary.^2) .* (E_secondary .< (E_primary - E_ionization) / 2)

            return secondary_e
        elseif s_or_c == "c"
            # Calculating the spectra of the degrading primary e-
            iLevel = findmin(abs.(E_ionization .- E_ionizations))[2]
            iPrimary = findmin(abs.(E4Q .- E_primary))[2]
            primary_e = Q[iPrimary, 1:length(E_secondary), iLevel]

            return primary_e
        end
    end
end


let Q = [], E4Q , E_ionizations
    global function cascading_O(E_secondary, E_primary, E_ionization, s_or_c)

        #=== WHERE DOES THAT COME FROM ? ===#
        E_parameters = [100, 200, 500, 1000, 2000]
        B_p = [7.18, 4.97, 2.75, 1.69, 1.02] .* 1e-22
        A_p = [12.6, 13.7, 14.1, 14.0, 13.7]
        if (minimum(E_parameters) < E_primary) & (E_primary < maximum(E_parameters))
            A = linear_interpolation(E_parameters, A_p)(E_primary)
            B = linear_interpolation(E_parameters, B_p)(E_primary)
        elseif E_primary <= minimum(E_parameters)
            A = A_p[1]
            B = B_p[1]
        else
            A = A_p[end]
            B = B_p[end]
        end
        A = A*1.25
        #===================================#

        # First we check if the cascading spectra has already been loaded by checking the
        # static/persistent local variables (namely Q, E4Q and E_ionizations)
        # check1, check2, check3 = is_cascading_loaded(Q, E4Q, E_secondary)

        # if check1 || check2 || check3
        if isempty(Q) || length(E_secondary) > length(E4Q) || E4Q[1:length(E_secondary)] != E_secondary
            # then we try to find a cascading spectra file with matching energy grid E_secondary
            found_them = 0
            cascading_files = readdir(pkgdir(AURORA, "internal_data", "e_cascading", "O"))
            for i1 in eachindex(cascading_files)
                if !isdir(cascading_files[i1])
                    try
                        filename = pkgdir(AURORA, "internal_data", "e_cascading", "O", cascading_files[i1])
                        file = matopen(filename)
                        E4Q = read(file, "E4Q")
                        close(file)
                        if length(E_secondary) <= length(E4Q) && E4Q[1:length(E_secondary)] == E_secondary
                            println("Loading cascading-matrices from file: ", cascading_files[i1])
                            file = matopen(filename)
                            Q = read(file, "Q")
                            E_ionizations = read(file, "E_ionizations")
                            close(file)
                            found_them = 1
                        break
                        end
                    catch
                    end
                end
            end
            if found_them == 0
                println("Could not find file with matching energy grid.")
                println("Starting to calculate the requested cascading-matrices.")

                # Precalculating the degrading spectra
                E_ionizations = [13.618, 16.9, 18.6, 28.5]
                Q = zeros(length(E_secondary), length(E_secondary), length(E_ionizations))
                E4Q = E_secondary
                dE = diff(E_secondary); dE = dE[[1:end; end]]
                println("Pre-calculating all energy-degradations for e - O-ionizing collisions.")
                println("Starting at ", Dates.format(now(), "HH:MM:SS"))
                for i1 in length(E_ionizations):-1:1
                    iLim = findfirst(x -> x > E_ionizations[i1], E_secondary)

                    # Lower boundary of primary electrons, This correctly accounts for the
                    # limit when Edeg[i3 + 1] is larger than Eprimary[i2 + 1] - Eionizations[i1],
                    # cutting a corner off the [Ei3 - (Ei3 + dEi3)] x [Ei2 + dEi2] square
                    # we're integrating over.
                    Emin(E) = max(Eprime[1], E + E_ionizations[i1])

                    # Change variable E[2] so that we can integrate 'fun' over an hypercube later
                    Eprime = Vector{Float64}    # we have to define an empty Eprime here so it can be
                                                # used in the function definition
                    var_hypercube(E1, E2) = Emin(E1) .+ E2 .* (Eprime[end] .- Emin(E1))
                    # Compute the jacobian associated to the change of variable above
                    jacobian(E1) = (Eprime[end] .- Emin(E1))

                    # Function fun... TODO
                    fun(E) = B ./
                            (1 + ((var_hypercube(E[1], E[2]) .- E_ionizations[i1] .- E[1]) / A).^(5/3)) .*
                            jacobian(E[1])

                    p = Progress(length(E_secondary) - iLim; desc=string("Doing level ", i1), color=:blue);
                    for i2 in iLim:length(E_secondary)

                        #=== WHERE DOES THAT COME FROM ? ===#
                        if (minimum(E_parameters) < E_secondary[i2]) & (E_secondary[i2] < maximum(E_parameters))
                            A = linear_interpolation(E_parameters, A_p)(E_secondary[i2])
                            B = linear_interpolation(E_parameters, B_p)(E_secondary[i2])
                        elseif E_secondary[i2] <= minimum(E_parameters)
                            A = A_p[1]
                            B = B_p[1]
                        else
                            A = A_p[end]
                            B = B_p[end]
                        end
                        #===================================#

                        iHalf = findlast(x -> x < (E_secondary[i2] - E_ionizations[i1]) / 2, E_secondary)
                        if !isnothing(iHalf) # otherwise we go to next i2
                            Eprime = [E_secondary[i2], E_secondary[i2] + dE[i2]]

                            for i3 in iHalf:(i2 - 1)
                                # Correct upper Edeg integration limit, this cuts the region into a
                                # triangle when needed.
                                Edegmax = min(E_secondary[i3] + dE[i3], E_secondary[i2] + dE[i2] - E_ionizations[i1])
                                # This if-condition just makes sure that we have the physically meaningful
                                # limits, i.e. we integrate in increasing Edeg - which is not otherwise
                                # guaranteed.
                                if Edegmax > E_secondary[i3]
                                    Q[i2, i3, i1] = hcubature(fun, [E_secondary[i3], 0],
                                                                    [Edegmax, 1], rtol = 1e-4)[1]
                                end
                            end
                        end
                        next!(p)
                    end
                    # println("Done with level ", i1, " at: ", Dates.format(now(), "HH:MM:SS"))
                end

                # Save the results for future use
                filename = pkgdir(AURORA, "internal_data", "e_cascading", "O",
                                    string("CascadingSpecOionization_",
                                    Dates.format(now(), "yyyymmdd-HHMMSS"),
                                    ".mat"))
                file = matopen(filename, "w")
                write(file, "Q", Q)
                write(file, "E4Q", E4Q)
                write(file, "E_ionizations", E_ionizations)
                close(file)
            end
        end



        if s_or_c == "s"
            # Calculating the spectra of the secondary e-
            secondary_e = B ./ (1 .+ (E_secondary ./ A).^(5/3)) .* (E_secondary .< (E_primary - E_ionization) / 2)

            return secondary_e
        elseif s_or_c == "c"
            # Calculating the spectra of the degrading primary e-
            iLevel = findmin(abs.(E_ionization .- E_ionizations))[2]
            iPrimary = findmin(abs.(E4Q .- E_primary))[2]
            primary_e = Q[iPrimary, 1:length(E_secondary), iLevel]

            return primary_e
        end
    end
end
