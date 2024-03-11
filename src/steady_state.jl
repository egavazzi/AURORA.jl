using KLU
using LinearAlgebra
using SparseArrays

function Steady_state(h_atm, μ, A, B, D, Q, Ie_top)
    Ie = Array{Float64}(undef, length(h_atm) * length(μ))

    # Spatial differentiation matrices, for up and down streams
    # Here we tuck on a fictious height one tick below lowest height and on tick above highest
    # height just to make it possible to use the diff function.
    h4diffu = [h_atm[1] - (h_atm[2] - h_atm[1]) ; h_atm]
    h4diffd = [h_atm ; h_atm[end] + (h_atm[end] - h_atm[end-1])]
    Ddz_Up   = spdiagm(-1 => -1 ./ (2 .* diff(h4diffu[2:end])),
                        0 =>  1 ./ (2 .* diff(h4diffu[1:end])))
    Ddz_Down = spdiagm( 0 => -1 ./ (2 .* diff(h4diffd[1:end])),
                        1 =>  1 ./ (2 .* diff(h4diffd[1:end-1])))

    # Diffusion operator
    Ddiffusion = d2M(h_atm)
    Ddiffusion[1, 1] = 0

    # Building the CN matrices
    Nz = length(h_atm)
    row_l = Vector{Int64}() # maybe using sizehint could help?
    col_l = Vector{Int64}()
    val_l = Vector{Float64}()
    for i1 in axes(B, 2)
        for i2 in axes(B, 2)
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
                    tmp_lhs =   μ[i1] .* Ddz_Down .+ Diagonal(A) .- D[i1] .* Ddiffusion .+ Diagonal(-B_tmp)
                    tmp_lhs[[1, end], :] .= 0
                    tmp_lhs[end, end] = 1
                else            # upward fluxes
                    tmp_lhs =   μ[i1] .* Ddz_Up .+ Diagonal(A) .- D[i1] .* Ddiffusion .+ Diagonal(-B_tmp)
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
    Mlhs = sparse(row_l, col_l, val_l)
    dropzeros!(Mlhs)    # for the performance of the next calculations

    index_top_bottom = sort(vcat(1:length(h_atm):(length(μ)*length(h_atm)),
                            length(h_atm):length(h_atm):(length(μ)*length(h_atm))))

    AAA = klu(Mlhs)

    I_top_bottom = (Ie_top * [0, 1]')'
    Q_local = copy(Q)
    Q_local[index_top_bottom] = I_top_bottom[:]
    Ie = AAA \ Q_local


    Ie[Ie .< 0] .= 0; # the fluxes should never be negative
    return Ie
end



using MAT
using ProgressMeter
using Dates
using Term
function steady_state_transport(altitude_max, θ_lims, E_max, B_angle_to_zenith,
    msis_file, iri_file, root_savedir, name_savedir, INPUT_OPTIONS)
    ## Get atmosphere
    println("Calling Matlab for the setup...")
    h_atm, ne, Te, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
    μ_scatterings = setup_new(altitude_max, θ_lims, E_max, msis_file, iri_file);

    ## Initialise
    I0 = zeros(length(h_atm) * length(μ_center), length(E));    # starting e- flux profile

    ## Load incoming flux
    if INPUT_OPTIONS.input_type == "from_old_matlab_file"
        Ie_top = Ie_top_from_old_matlab_file(1:1:1, E, 1, μ_center, INPUT_OPTIONS.input_file);
    elseif INPUT_OPTIONS.input_type == "from_file"
        Ie_top = Ie_top_from_file(1:1:1, E, μ_center, 1, INPUT_OPTIONS.input_file)
    elseif INPUT_OPTIONS.input_type == "flickering"
        Ie_top = Ie_top_flickering(1:1:1, E, dE, 1, μ_center, h_atm,
                                    μ_scatterings.BeamWeight, INPUT_OPTIONS.IeE_tot,
                                    INPUT_OPTIONS.z₀, INPUT_OPTIONS.E_min, INPUT_OPTIONS.f,
                                    INPUT_OPTIONS.Beams, INPUT_OPTIONS.modulation)
    elseif INPUT_OPTIONS.input_type == "constant_onset"
        Ie_top = Ie_top_constant(1:1:1, E, dE, 1, μ_center, h_atm,
                                μ_scatterings.BeamWeight, INPUT_OPTIONS.IeE_tot,
                                INPUT_OPTIONS.z₀, INPUT_OPTIONS.E_min, INPUT_OPTIONS.Beams,
                                INPUT_OPTIONS.t0, INPUT_OPTIONS.t1)
    end

    ## Calculate the phase functions and put them in a Tuple
    phaseN2e, phaseN2i = phase_fcn_N2(μ_scatterings.theta1, E);
    phaseO2e, phaseO2i = phase_fcn_O2(μ_scatterings.theta1, E);
    phaseOe, phaseOi = phase_fcn_O(μ_scatterings.theta1, E);
    phase_fcn_neutrals = ((phaseN2e, phaseN2i), (phaseO2e, phaseO2i), (phaseOe, phaseOi));
    cascading_neutrals = (cascading_N2, cascading_O2, cascading_O) # tuple of functions


    ## Create the folder to save the data to
    if isempty(root_savedir) || !occursin(r"[^ ]", root_savedir)
        # if root_savedir is empty or contains only "space" characters, we use "backup/" as a name
        root_savedir = "backup"
    end
    if isempty(name_savedir) || !occursin(r"[^ ]", name_savedir)
        # if name_savedir is empty or contains only "space" characters, we use the current
        # date and time as a name
        name_savedir = string(Dates.format(now(), "yyyymmdd-HHMM"))
    end
    savedir = pkgdir(AURORA, "data", root_savedir, name_savedir)
    if isdir(savedir) && (filter(startswith("IeFlickering-"), readdir(savedir)) |> length) > 0
        # throw a warning if name_savedir exists and if it already contains results
        print("\n", @bold @red "WARNING!")
        print(@bold " '$savedir' ")
        println(@bold @red "already exists, any results stored in it will be overwritten.")
        # println(@bold @red "already exists, the experiment is aborted.")
        # return
    else
        if ~isdir(pkgdir(AURORA, "data", root_savedir)) # check if the root_savedir exists
            mkdir(pkgdir(AURORA, "data", root_savedir)) # if not, creates it
        end
        mkpath(savedir)
    end
    print("\n", @bold "Results will be saved at $savedir \n")

    ## And save the simulation parameters in it
    save_parameters(altitude_max, θ_lims, E_max, B_angle_to_zenith, 1:1:1, 1:1:1, 1,
        0, INPUT_OPTIONS, savedir)
    save_neutrals(h_atm, n_neutrals, ne, Te, savedir)

    # Initialize arrays for the ionization collisions part of the energy degradation
    Ionization_matrix = [zeros(length(h_atm) * length(μ_center), 1) for _ in 1:15]
    Ionizing_matrix = [zeros(length(h_atm) * length(μ_center), 1) for _ in 1:15]
    secondary_vector = [zeros(length(E)) for _ in 1:15]
    primary_vector = [zeros(length(E)) for _ in 1:15]


    ## No n_loop here
    Q  = zeros(length(h_atm) * length(μ_center), 1, length(E));
    Ie = zeros(length(h_atm) * length(μ_center), 1, length(E));

    D = make_D(E, dE, θ_lims);
    # Extract the top flux for the current loop
    Ie_top_local = Ie_top[:, 1, :];

    p = Progress(length(E), desc=string("Calculating flux"))
    # Looping over energy
    for iE in length(E):-1:1
        A = make_A(n_neutrals, σ_neutrals, ne, Te, E, dE, iE);

        B, B2B_inelastic_neutrals = make_B(n_neutrals, σ_neutrals, E_levels_neutrals,
                                phase_fcn_neutrals, dE, iE, μ_scatterings.Pmu2mup,
                                μ_scatterings.BeamWeight_relative, μ_scatterings.theta1);

        # Compute the flux of e-
        Ie[:, 1, iE] = Steady_state(h_atm ./ cosd(B_angle_to_zenith), μ_center, A, B,
                                    D[iE, :],  Q[:, 1, iE], Ie_top_local[:, iE])

        # Update the cascading of e-
        update_Q!(Q, Ie, h_atm, 1:1:1, ne, Te, n_neutrals, σ_neutrals, E_levels_neutrals, B2B_inelastic_neutrals,
                    cascading_neutrals, E, dE, iE, μ_scatterings.BeamWeight, μ_center,
                    Ionization_matrix, Ionizing_matrix, secondary_vector, primary_vector)

        next!(p)
    end

    save_results(Ie, E, 1:1:1, μ_lims, h_atm, I0, μ_scatterings, 1, 1, savedir)

end
