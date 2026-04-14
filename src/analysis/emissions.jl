using MAT: matopen, matread
using ProgressMeter: Progress, next!


# ======================================================================================== #
#                                  OPTICAL EMISSIONS                                     #
# ======================================================================================== #

"""
    make_volume_excitation_file(directory_to_process)

Reads into a folder `directory_to_process` containing results from an AURORA.jl simulation,
loads the particle flux `Ie` (#e⁻/m²/s), and calculates the volume-excitation-rates. For
prompt emissions, volume-excitation-rates correspond also to volume-emission-rates.

# Calling
`make_volume_excitation_file(directory_to_process)`

# Inputs
- `directory_to_process`: absolute or relative path to the simulation directory to process.
"""
function make_volume_excitation_file(directory_to_process)
    ## Find the files to process
    files = readdir(directory_to_process, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]
    # The files are sorted in lexicographical order, so IeFlickering-100.mat will be loaded
    # before "IeFlickering-11.mat. We fix that with the following line which sorts them by
    # the number in the filename.
    sort!(files_to_process, by = x -> parse(Int, match(r"IeFlickering-(\d+)\.mat", basename(x))[1]))

    if isempty(files_to_process)
        @warn "No simulation results found in $directory_to_process. Skipping volume excitation calculations."
        return nothing
    end

    ## Load simulation grid
    f = matopen(files_to_process[1])
        z = read(f, "h_atm")
        E_centers = read(f, "E_centers")
    close(f)

    ## Load simulation neutral densities
    data = matread(joinpath(directory_to_process, "neutral_atm.mat"))
    nO = data["nO"]
    nO2 = data["nO2"]
    nN2 = data["nN2"]

    ## Load the emission cross-sections
    σ_4278 = excitation_4278(E_centers)
    σ_6730 = excitation_6730_N2(E_centers)
    σ_7774_O = excitation_7774_O(E_centers)
    σ_7774_O2 = excitation_7774_O2(E_centers)
    σ_8446_O = excitation_8446_O(E_centers)
    σ_8446_O2 = excitation_8446_O2(E_centers)
    σ_O1D = excitation_O1D(E_centers)
    σ_O1S = excitation_O1S(E_centers)
    ## Load/calculate the ionization cross-sections
    σ_N2, σ_O2, σ_O = load_cross_sections(E_centers) # load the cross-sections
    N2_levels, O2_levels, O_levels = load_excitation_threshold() # load the energy levels
    σ_Oi = σ_O' * O_levels[:, 2] # basically do sum(cross_section_for_each_reaction * number_of_ionizations_per_reaction)
    σ_O2i = σ_O2' * O2_levels[:, 2]
    σ_N2i = σ_N2' * N2_levels[:, 2]

    ## Initialize arrays to store the results for each time-slice
    Q4278 = Vector{Matrix{Float64}}()
    Q6730 = Vector{Matrix{Float64}}()
    Q7774_O = Vector{Matrix{Float64}}()
    Q7774_O2 = Vector{Matrix{Float64}}()
    Q7774 = Vector{Matrix{Float64}}()
    Q8446_O = Vector{Matrix{Float64}}()
    Q8446_O2 = Vector{Matrix{Float64}}()
    Q8446 = Vector{Matrix{Float64}}()
    QO1D = Vector{Matrix{Float64}}()
    QO1S = Vector{Matrix{Float64}}()
    QOi = Vector{Matrix{Float64}}()
    QO2i = Vector{Matrix{Float64}}()
    QN2i = Vector{Matrix{Float64}}()
    t = Float64[]

    n_files = length(files_to_process)
    p = Progress(n_files; desc=string("Processing data"), dt=1.0, color=:blue)
    ## Loop over the files
    for (i_file, file) in enumerate(files_to_process)
        ## Load simulation results for current file.
        if i_file == 1
            f = matopen(file)
                Ie_ztE = read(f, "Ie_ztE")
                t_local = read(f, "t_run")
            close(f)
        else
            f = matopen(file)
            @views Ie_ztE = read(f, "Ie_ztE")[:, 2:end, :]
            t_local = read(f, "t_run")[2:end]
            close(f)
        end

        ## Sum Ie over the beams
        #=
        I had the idea of initializing Ie_ztE_omni outside of the file-reading loop and then only
        fill it with zeros at the beginning of each new loop, to minimize allocations and GC.
        But then I realized that the current version of the code with the initialization inside
        the loop allows for the Ie_ztE slices to have different n_t sizes.
        Currently (v0.4.2), the Ie_ztE saved in each IeFlickering-XX.mat file all have the
        same size so it is not of big use. But this is not a bottleneck anyway and it allows
        for using time slices with different n_t sizes in the future. // EG 20241221
        =#
        n_z = length(z)
        n_μ = size(Ie_ztE, 1) ÷ n_z # (÷ returns an Int)
        n_t = size(Ie_ztE, 2)
        n_E = size(Ie_ztE, 3)
        Ie_ztE_omni = zeros(n_z, n_t, n_E)
        @views for i_μ in 1:n_μ
            Ie_ztE_omni .+= Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :] # bottleneck, ~75% of time spent here
        end

        # This was an attempt to make things faster with @turbo, but I didn't observe big
        # differences when benchmarking on my machine.
        # @turbo for i_E in 1:n_E
        #     for i_t in 1:n_t
        #         for i_μ in 1:n_μ
        #             for i_z in 1:n_z
        #                 Ie_ztE_omni[i_z, i_t, i_E] += Ie_ztE[(i_μ - 1) * n_z + i_z, i_t, i_E]
        #             end
        #         end
        #     end
        # end

        ## Calculate Q (volume-excitation-rate) for various optical emissions
        Q4278_local = calculate_volume_excitation(z, t_local, Ie_ztE_omni, σ_4278, nN2)
        Q6730_local = calculate_volume_excitation(z, t_local, Ie_ztE_omni, σ_6730, nN2)
        Q7774_O_local = calculate_volume_excitation(z, t_local, Ie_ztE_omni, σ_7774_O, nO)
        Q7774_O2_local = calculate_volume_excitation(z, t_local, Ie_ztE_omni, σ_7774_O2, nO2)
        Q7774_local = Q7774_O_local + Q7774_O2_local
        Q8446_O_local = calculate_volume_excitation(z, t_local, Ie_ztE_omni, σ_8446_O, nO)
        Q8446_O2_local = calculate_volume_excitation(z, t_local, Ie_ztE_omni, σ_8446_O2, nO2)
        Q8446_local = Q8446_O_local + Q8446_O2_local
        QO1D_local = calculate_volume_excitation(z, t_local, Ie_ztE_omni, σ_O1D, nO) # quenching is not taken into account
        QO1S_local = calculate_volume_excitation(z, t_local, Ie_ztE_omni, σ_O1S, nO) # quenching is not taken into account
        # Calculate Q (volume-excitation-rate) for ionizations
        QOi_local = calculate_volume_excitation(z, t_local, Ie_ztE_omni, σ_Oi, nO)
        QO2i_local = calculate_volume_excitation(z, t_local, Ie_ztE_omni, σ_O2i, nO2)
        QN2i_local = calculate_volume_excitation(z, t_local, Ie_ztE_omni, σ_N2i, nN2)

        ## Push the newly calculated Q_local for the current time-slice into a vector
        # We get something like Q4278 = [[n_z, n_t1], [n_z, n_t2], ...]
        push!(Q4278, Q4278_local)
        push!(Q6730, Q6730_local)
        push!(Q7774_O, Q7774_O_local)
        push!(Q7774_O2, Q7774_O2_local)
        push!(Q7774, Q7774_local)
        push!(Q8446_O, Q8446_O_local)
        push!(Q8446_O2, Q8446_O2_local)
        push!(Q8446, Q8446_local)
        push!(QO1D, QO1D_local)
        push!(QO1S, QO1S_local)
        push!(QOi, QOi_local)
        push!(QO2i, QO2i_local)
        push!(QN2i, QN2i_local)
        append!(t, t_local)

        next!(p)
    end

    ## Concatenate along time
    # We get Q4278 = [n_z, n_t]
    Q4278 = reduce(hcat, Q4278)
    Q6730 = reduce(hcat, Q6730)
    Q7774_O = reduce(hcat, Q7774_O)
    Q7774_O2 = reduce(hcat, Q7774_O2)
    Q7774 = reduce(hcat, Q7774)
    Q8446_O = reduce(hcat, Q8446_O)
    Q8446_O2 = reduce(hcat, Q8446_O2)
    Q8446 = reduce(hcat, Q8446)
    QO1D = reduce(hcat, QO1D)
    QO1S = reduce(hcat, QO1S)
    QOi = reduce(hcat, QOi)
    QO2i = reduce(hcat, QO2i)
    QN2i = reduce(hcat, QN2i)

    ## Save results
    savefile = joinpath(directory_to_process, "Qzt_all_L.mat")
    f = matopen(savefile, "w")
        write(f, "h_atm", z)
        write(f, "t", t)
        write(f, "Q4278", Q4278)
        write(f, "Q6730", Q6730)
        write(f, "Q7774_O", Q7774_O)
        write(f, "Q7774_O2", Q7774_O2)
        write(f, "Q7774", Q7774)
        write(f, "Q8446_O", Q8446_O)
        write(f, "Q8446_O2", Q8446_O2)
        write(f, "Q8446", Q8446)
        write(f, "QO1D", QO1D)
        write(f, "QO1S", QO1S)
        write(f, "QOi", QOi)
        write(f, "QO2i", QO2i)
        write(f, "QN2i", QN2i)
    close(f)

    println("Volume excitation rates saved in $savefile")

    return nothing
end


"""
    calculate_volume_excitation(z, t, Ie_ztE_omni, σ, n)

Calculate the volume-excitation-rate for an excitation of interest, produced by the electron
flux `Ie_ztE_omni` that is summed over the beams (omnidirectional).

The excitation of interest is chosen through the cross-section `σ` given to the function.
Note that the neutral density `n` should match the excitation of interest (e.g. use nN2 when
calculating the volume-excitation-rate of the 4278Å optical emission).

# Calling
`Q = calculate_volume_excitation(z, t, Ie_ztE_omni, σ, n)`

# Inputs
- `z`: altitude (m). Vector [n\\_z]
- `t`: time (s). Vector [n\\_t]
- `Ie_ztE_omni`: omnidirectional electron flux (#e⁻/m²/s). 3D array [n\\_z, n\\_t, n\\_E]
- `σ`: excitation cross-section (m⁻²). Vector [n\\_E]
- `n`: density of exciteable atmospheric specie (m⁻³). Vector [n\\_z]
"""
function calculate_volume_excitation(z, t, Ie_ztE_omni, σ, n)
    # Initialize
    n_z = length(z)
    n_t = length(t)
    Q = zeros(n_z, n_t)
    # Calculate Q for each time step
    @views for i_t in eachindex(t)
        Q[:, i_t] .= (Ie_ztE_omni[:, i_t, :] * σ) .* n
    end

    return Q
end


"""
    make_column_excitation_file(directory_to_process)

Reads into a folder `directory_to_process` containing results from an AURORA.jl simulation,
loads the volume excitation rates `Q_XXXX` (#excitation/m³/s) contained in the file `Qzt_all_L.mat`
and integrate them in height, taking into account the finite speed of light.

The calculated colum-integrated excitation rates are saved to a file named `I_lambda_of_t.mat`.
The column-integrated excitation rates are named "I\\_4278, I\\_6730, ...". They are all vectors
in time (length n\\_t), and have units of (#excitation/m²/s).

Note that the function `make_volume_excitation_file()` needs to be run before this one, as
we need the file `Qzt_all_L.mat` with the volume excitation rates.

# Calling
`make_column_excitation_file(directory_to_process)`

# Inputs
- `directory_to_process`: absolute or relative path to the simulation directory to process.
"""
function make_column_excitation_file(directory_to_process)
    ## Find the files to process
    Q_file = joinpath(directory_to_process, "Qzt_all_L.mat")

    ## Load volume-excitation-rates
    data = matread(Q_file)
    z = data["h_atm"]
    t = data["t"]
    Q4278 = data["Q4278"]
    Q6730 = data["Q6730"]
    Q7774 = data["Q7774"]
    Q7774_O = data["Q7774_O"]
    Q7774_O2 = data["Q7774_O2"]
    Q8446 = data["Q8446"]
    Q8446_O = data["Q8446_O"]
    Q8446_O2 = data["Q8446_O2"]
    QO1D = data["QO1D"]
    QO1S = data["QO1S"]

    ## Integrate in altitude, taking into account finite photon velocity and path length.
    I_4278 = q2colem(t, z, Q4278)
    I_6730 = q2colem(t, z, Q6730)
    I_7774 = q2colem(t, z, Q7774)
    I_7774_O = q2colem(t, z, Q7774_O)
    I_7774_O2 = q2colem(t, z, Q7774_O2)
    I_8446 = q2colem(t, z, Q8446)
    I_8446_O = q2colem(t, z, Q8446_O)
    I_8446_O2 = q2colem(t, z, Q8446_O2)
    I_O1D = q2colem(t, z, QO1D)
    I_O1S = q2colem(t, z, QO1S)

    ## Save results
    savefile = joinpath(directory_to_process, "I_lambda_of_t.mat")
    f = matopen(savefile, "w")
        write(f, "t", t)
        write(f, "I_4278", I_4278)
        write(f, "I_6730", I_6730)
        write(f, "I_7774", I_7774)
        write(f, "I_7774_O", I_7774_O)
        write(f, "I_7774_O2", I_7774_O2)
        write(f, "I_8446", I_8446)
        write(f, "I_8446_O", I_8446_O)
        write(f, "I_8446_O2", I_8446_O2)
        write(f, "I_O1D", I_O1D)
        write(f, "I_O1S", I_O1S)
    close(f)

    println("Column excitation rates saved in $savefile")

    return nothing
end

using Integrals: SampledIntegralProblem, TrapezoidalRule, solve
using Interpolations: interpolate, extrapolate, Gridded, Linear
"""
    q2colem(t::Vector, z, Q, A = 1, τ = ones(length(z)))

Integrate the volume-excitation-rate (#exc/m³/s) to column-excitation-rate (#exc/m²/s).

Takes into account the time-delay between light emitted at different altitudes. Photons
emitted at at altitude of 200km will arrive at the detector 100e3/3e8 = 0.333 ms later than
electrons emitted at an altitude of 100km. This is a small time-shift, but it is close to
the time-differences corresponding to the phase-shifts between auroral emissions varying at
~10Hz.

The einstein coefficient `A` and effective lifetime `τ` are optional (equal to one by default).

# Calling
`I = q2colem(t, z, Q, A, τ)`

# Inputs
- `z`: altitude (m). Vector [n\\_z]
- `t`: time (s). Vector [n\\_t]
- `Q`: volume-excitation-rate (#exc/m³/s) of the wavelength of interest. 2D array [n\\_z, n\\_t]
- `A`: einstein coefficient (s⁻¹). Scalar (Float or Int)
- `τ`: effective lifetime (s). Vector [n\\_z].

# Output
- `I`: integrated column-excitation-rate (#exc/m²/s) of the wavelength of interest. Vector [n\\_t]
"""
function q2colem(t::Vector, z, Q, A = 1, τ = ones(length(z)))

    ## Define constant
    c = 2.99792458e8 # speed of light (m/s)

    ## Apply the effective lifetime and the Einstein coefficient
    Q = Q .* τ .* A

    ## Create the 2D interpolator
    # To explain how this function is working, let's follow an example.
    # Suppose we have the following Q matrix:
    #     julia> Q
    #     5×5 Matrix{Float64}:
    #     0.65745   0.0652896  0.313073   0.4075     0.811552
    #     0.780153  0.530831   0.546205   0.575431   0.707561
    #     0.573467  0.555628   0.567592   0.88814    0.728269
    #     0.865137  0.141949   0.0853166  0.0506817  0.65735
    #     0.707699  0.646255   0.959993   0.932128   0.956778
    # given over the following time and height vectors:
    #     julia> t
    #     5-element Vector{Int64}:
    #         1
    #         2
    #         3
    #         4
    #         5
    #     julia> z
    #     5-element Vector{Int64}:
    #         1
    #         2
    #         3
    #         4
    #         5
    nodes = (z, t)
    itp = interpolate(nodes, Q, Gridded(Linear()))
    itp = extrapolate(itp, 0.0) # allows for extrapolation but set those values to 0.0

    ## Shift data in time
    # We shift the values in time.
    # For each height (index i) and time (index j) position in the Q matrix, we calculate at
    # what time the corresponding photons will reach the bottom of the height column. That new
    # time is given by the formula `t[j] - (z[i] - z[1]) / c`. We use that new time as
    # input to the interpolator `itp`,  which has for effect to "shift" that value in time.
    # Continuing our example from above, and taking c = 1 for simplicity, we get
    #     julia> I = [itp(z[i], (t[j] - (z[i] - z[1]) / 1)) for i in eachindex(z), j in eachindex(t)]
    #     5×5 Matrix{Float64}:
    #     0.65745  0.0652896  0.313073  0.4075    0.811552
    #     0.0      0.780153   0.530831  0.546205  0.575431
    #     0.0      0.0        0.573467  0.555628  0.567592
    #     0.0      0.0        0.0       0.865137  0.141949
    #     0.0      0.0        0.0       0.0       0.707699
    # We can see that the values have been shifted/delayed in time.
    # For the first time t = 1, only photons from the height z = 1 are arriving. For the
    # time t = 2, photons from the height z = 2 that were emitted at time t = 1 are now
    # arriving. Etc. etc.

    # We used a 2D interpolation. As the height are unchanged, we could have also done a 1D
    # interpolation for each height. We use the 2D interpolation for convenience as it allows
    # us to have only one interpolatior object (instead of one per height). This is also
    # closest to the method used in the legacy Matlab code from which this function is
    # inspired.
    I = [itp(z[i], (t[j] - (z[i] - z[1]) / c)) for i in eachindex(z), j in eachindex(t)]

    ## Integrate in height
    # Now that we have the matrix I which contains the information about the number of
    # photons arriving at the bottom of the column for each height and time steps, we can
    # integrate it along the height.
    problem = SampledIntegralProblem(I, z; dim=1)
    method = TrapezoidalRule()
    I_lambda = solve(problem, method)

    return I_lambda.u
end

## Steady-state version
"""
    q2colem(t::Real, z, Q, A = 1, τ = ones(length(z)))

Same as above, except time is now a scalar (steady-state results). This is just a simple
integration in height.
"""
function q2colem(t::Real, z, Q, A = 1, τ = ones(length(z)))
    ## Apply the effective lifetime and the Einstein coefficient
    Q = Q .* τ .* A

    ## Simple 1D integration
    problem = SampledIntegralProblem(Q, z; dim=1)
    method = TrapezoidalRule()
    I_lambda = solve(problem, method)

    return I_lambda.u
end


# ======================================================================================== #
#                         AuroraSimulation convenience wrappers                          #
# ======================================================================================== #

"""
    make_volume_excitation_file(sim::AuroraSimulation)

Convenience wrapper that calls [`make_volume_excitation_file`](@ref) on `sim.savedir`.
"""
make_volume_excitation_file(sim::AuroraSimulation) = make_volume_excitation_file(sim.savedir)

"""
    make_column_excitation_file(sim::AuroraSimulation)

Convenience wrapper that calls [`make_column_excitation_file`](@ref) on `sim.savedir`.
"""
make_column_excitation_file(sim::AuroraSimulation) = make_column_excitation_file(sim.savedir)
