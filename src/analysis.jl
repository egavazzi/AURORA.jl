using MAT: matopen, matread
using ProgressMeter: Progress, next!


"""
    make_psd_file(directory_to_process; compute=:both, vpar_edges=nothing, output_prefix="psd")

Reads into a folder `directory_to_process` containing results from an AURORA.jl simulation,
loads particle flux files `IeFlickering-NN.mat`, converts them into phase-space density, and
writes one output file per input as `psd/psd-NN.mat`.

# Calling
`make_psd_file(directory_to_process; compute=:both, vpar_edges=nothing, output_prefix="psd")`

# Inputs
- `directory_to_process`: absolute or relative path to the simulation directory to process.

# Keyword Arguments
- `compute`: one of `:f_only`, `:F_only`, or `:both`.
- `vpar_edges`: custom `v_parallel` bin edges [m/s] used when computing `F`.
    If `nothing`, an automatic symmetric uniform interval grid is used, spanning
    `[-maximum(v), maximum(v)]` with an edge at `v_parallel = 0`.
- `output_prefix`: output filename prefix used as `<output_prefix>-NN.mat`.
"""
function make_psd_file(
    directory_to_process;
    compute::Symbol = :both,
    vpar_edges::Union{Nothing, AbstractVector} = nothing,
    output_prefix::AbstractString = "psd",
)
    println("Converting Ie to PSD.")
    files = readdir(directory_to_process, join = true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]
    sort!(
        files_to_process,
        by = x -> parse(Int, match(r"IeFlickering-(\d+)\.mat", basename(x))[1]),
    )

    if isempty(files_to_process)
        @warn "No IeFlickering files found in $directory_to_process."
        return nothing
    end

    n_files = length(files_to_process)
    psd_dir = joinpath(directory_to_process, "psd")
    mkpath(psd_dir)

    for (i_file, file) in enumerate(files_to_process)
        match_result = match(r"IeFlickering-(\d+)\.mat", basename(file))
        @assert !isnothing(match_result) "Unexpected input filename: $file"
        file_id = match_result[1]

        print("\rProcessing PSD: $(i_file)/$(n_files)")
        flush(stdout)

        result = make_psd_from_AURORA(file; compute = compute, vpar_edges = vpar_edges)
        outfile = joinpath(psd_dir, "$(output_prefix)-$(file_id).mat")
        write_psd_result(outfile, result)
    end

    println()

    return nothing
end


function write_psd_result(outfile::AbstractString, res)
    matopen(outfile, "w") do io
        if hasproperty(res, :f)
            write(io, "f", res.f)
        end

        if hasproperty(res, :F)
            write(io, "F", res.F)
            write(io, "vpar_edges", res.vpar_edges)
            write(io, "vpar_centers", res.vpar_centers)
            write(io, "dvpar", res.Δvpar)
        end

        write(io, "v", res.v)
        write(io, "v_par", res.v_par)
        write(io, "v_perp", res.v_perp)
        write(io, "dE_J", res.ΔE_J)
        write(io, "BeamWeight", res.BeamWeight)
        write(io, "mu_center", res.μ_center)
        write(io, "E", res.E)
        write(io, "t_run", res.t_run)
        write(io, "h_atm", res.h_atm)
        write(io, "mu_lims", res.μ_lims)
    end
end










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
    t = Vector{}()

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
        push!(t, t_local)

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
    t = reduce(vcat, t)

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
`Q = calculate_volume_excitation(z, t, Ie_ztE, σ, n)``

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
`make_current_file(directory_to_process)`

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
#                                  OTHER FUNCTIONS                                     #
# ======================================================================================== #

"""
    make_Ie_top_file(directory_to_process)

Reads into a folder `directory_to_process` containing results from an AURORA.jl simulation
and extracts the particle flux `Ie` (#e⁻/m²/s) at the top of the ionosphere (i.e. at the
max altitude used in the simulation).

# Calling
`make_Ie_top_file(directory_to_process)`

# Inputs
- `directory_to_process`: absolute or relative path to the simulation directory to process.
"""
function make_Ie_top_file(directory_to_process)
    ## Find the files to process
    files = readdir(directory_to_process, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]
    # The files are sorted in lexicographical order, so IeFlickering-100.mat will be loaded
    # before "IeFlickering-11.mat. We fix that with the following line which sorts them by
    # the number in the filename.
    sort!(files_to_process, by = x -> parse(Int, match(r"IeFlickering-(\d+)\.mat", basename(x))[1]))

    ## Load simulation grid
    f = matopen(files_to_process[1])
        z = read(f, "h_atm")
        E_centers = read(f, "E_centers")
        ΔE = read(f, "dE")
        scattering_data = read(f, "mu_scatterings")
    close(f)
    n_z = length(z)
    Ω_beam = scattering_data["BeamWeight"]

    ## Initialize arrays to store the results for each time-slice
    Ie_top = Vector{Array{Float64, 3}}()
    t = Vector{}()

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

        ## Extract Ie at the top for each beam
        Ie_top_local = Ie_ztE[n_z:n_z:end, :, :]

        ## Push the Ie_top of the current time-slice into a vector
        # We get Ie_top = [[n_μ, n_t1, n_E], [n_μ, n_t2, n_E], ...]
        push!(Ie_top, Ie_top_local)
        push!(t, t_local)

        next!(p)
    end

    ## Concatenate along time
    # We get Ie_top = [n_μ, n_t, n_E]
    Ie_top = reduce(hcat, Ie_top)
    t = reduce(vcat, t)

    ## Play with the units
    Ie_top_raw = copy(Ie_top) # in #e-/m²/s
    Ie_top = Ie_top ./ reshape(ΔE, (1, 1, :)) ./ Ω_beam # in #e-/m²/s/eV/ster

    ## Save results
    savefile = joinpath(directory_to_process, "Ie_top.mat")
    f = matopen(savefile, "w")
        write(f, "E_centers", E_centers)
        write(f, "dE", ΔE)
        write(f, "BeamW", Ω_beam)
        write(f, "t", t)
        write(f, "Ie_top_raw", Ie_top_raw)
        write(f, "Ie_top", Ie_top)
    close(f)

    println("Top flux saved in $savefile")

    return nothing
end


"""
    make_current_file(directory_to_process)

Reads into a folder `directory_to_process` containing results from an AURORA.jl simulation,
loads the particle flux `Ie` (#e⁻/m²/s) and calculates the field-aligned current-density
and field-aligned energy-flux for each height and through time.

The following variables are saved to a file named *J.mat*:
- `J_up`: Field-aligned current-density in the upward direction. 2D array [n\\_z, n\\_t]
- `J_down`: Field-aligned current-density in the downward direction. 2D array [n\\_z, n\\_t]
- `E_up`: Field-aligned energy-flux (eV/m²/s) in the upward direction. 2D array [n\\_z, n\\_t]
- `E_down`: Field-aligned energy-flux (eV/m²/s) in the downward direction. 2D array [n\\_z, n\\_t]

# Calling
`make_current_file(directory_to_process)`

# Inputs
- `directory_to_process`: absolute or relative path to the simulation directory to process.
"""
function make_current_file(directory_to_process)
    ## Find the files to process
    files = readdir(directory_to_process, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]
    # The files are sorted in lexicographical order, so IeFlickering-100.mat will be loaded
    # before "IeFlickering-11.mat. We fix that with the following line which sorts them by
    # the number in the filename.
    sort!(files_to_process, by = x -> parse(Int, match(r"IeFlickering-(\d+)\.mat", basename(x))[1]))

    ## Load simulation grid
    f = matopen(files_to_process[1])
        z = read(f, "h_atm")
        E_centers = read(f, "E_centers")
        μ_lims = read(f, "mu_lims")
    close(f)
    μ_center = mu_avg(acosd.(μ_lims))

    ## Initialize arrays to store the results for each time-slice
    J_up = Vector{Array{Float64, 2}}()
    J_down = Vector{Array{Float64, 2}}()
    IeE_up = Vector{Array{Float64, 2}}()
    IeE_down = Vector{Array{Float64, 2}}()
    t = Vector{}()

    ## Define constant
    q_e = 1.602176620898e-19 # elementary charge (C)

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

        ## Calculate the field aligned currents and energy flux
        n_z = length(z)
        n_μ = length(μ_center)
        n_t = size(Ie_ztE, 2)
        J_up_local = zeros(n_z, n_t)
        J_down_local = zeros(n_z, n_t)
        IeE_up_local = zeros(n_z, n_t)
        IeE_down_local = zeros(n_z, n_t)
        @views for i_μ in 1:n_μ
            if μ_center[i_μ] > 0
                J_up_local .+= q_e * abs(μ_center[i_μ]) .* sum(Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :], dims=3)
                IeE_up_local .+= abs(μ_center[i_μ]) .* sum(Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :] .* reshape(E_centers, (1, 1, :)), dims=3)
            else
                J_down_local .+= q_e * abs(μ_center[i_μ]) .* sum(Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :], dims=3)
                IeE_down_local .+= abs(μ_center[i_μ]) .* sum(Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :] .* reshape(E_centers, (1, 1, :)), dims=3)
            end
        end

        ## Push the J of the current time-slice into a vector
        # We get J_up = [[n_z, n_t1], [n_z, n_t2], ...]
        push!(J_up, J_up_local)
        push!(J_down, J_down_local)
        push!(IeE_up, IeE_up_local)
        push!(IeE_down, IeE_down_local)
        push!(t, t_local)

        next!(p)
    end

    ## Concatenate along time
    # We get J_up = [n_z, n_t]
    J_up = reduce(hcat, J_up)
    J_down = reduce(hcat, J_down)
    IeE_up = reduce(hcat, IeE_up)
    IeE_down = reduce(hcat, IeE_down)
    t = reduce(vcat, t)

    ## Save results
    savefile = joinpath(directory_to_process, "J.mat")
    f = matopen(savefile, "w")
        write(f, "h_atm", z)
        write(f, "t", t)
        write(f, "J_up", J_up)
        write(f, "J_down", J_down)
        write(f, "IeE_up", IeE_up)
        write(f, "IeE_down", IeE_down)
    close(f)

    println("Currents saved in $savefile")

    return nothing
end


"""
    downsampling_fluxes(directory_to_process, downsampling_factor)

This function extracts `Ie` from the simulation results in `directory_to_process` and
downsample it in time.
For example: if `Ie` is given with a time step of 1ms and we use a `downsampling_factor` of
10, this function will extract the values of `Ie` with a time step of 10ms. It will then
save the results in a new subfolder called`downsampled_10x`, inside the `directory_to_process`.

# Calling
`downsampling_fluxes(directory_to_process, downsampling_factor)`

# Inputs
- `directory_to_process`: absolute or relative path to the simulation directory to process.
- `downsampling_factor`: downsampling factor for the time

# Outputs
The downsampled electron fluxes `Ie` will be saved in a subfolder inside the `directory_to_process`.
"""
function downsampling_fluxes(directory_to_process, downsampling_factor)
    files = readdir(directory_to_process, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]
    # The files are sorted in lexicographical order, so IeFlickering-100.mat will be loaded
    # before "IeFlickering-11.mat. We fix that with the following line which sorts them by
    # the number in the filename.
    sort!(files_to_process, by = x -> parse(Int, match(r"IeFlickering-(\d+)\.mat", basename(x))[1]))

    for j in files_to_process
        f = matopen(j)
            Ie = read(f, "Ie_ztE") # [n_μ * nz, nt, nE]
            E_centers = read(f, "E_centers")
            ΔE = read(f, "dE")
            t_run = read(f, "t_run")
            μ_lims = read(f, "mu_lims")
            z = read(f, "h_atm")
            I0 = read(f, "I0")
            scattering_data = read(f, "mu_scatterings")
        close(f)

        # downsample the data
        dt = diff(t_run)[1]
        println("The time-step from simulation is ", dt, "s.")
        new_dt = dt * downsampling_factor
        println("The time-step of the new file will be ", new_dt, "s.")
        Ie = Ie[:, 1:downsampling_factor:end, :]
        t_run = t_run[1:downsampling_factor:end]

        # create new subdir if it doesn't exist
        new_subdir = "downsampled_" * string(downsampling_factor) * "x"
        full_path_to_new_subdir = joinpath(directory_to_process, new_subdir)
        mkpath(full_path_to_new_subdir)

        # create new file
        new_filename = splitdir(j)[2][1:end-4] * "d.mat" # add a 'd' after the number to indicate that it is downsampled
        full_path_to_new_filename = joinpath(full_path_to_new_subdir, new_filename)
        file = matopen(full_path_to_new_filename, "w")
            write(file, "Ie_ztE", Ie)
            write(file, "E_centers", E_centers)
            write(file, "E_edges", E_edges)
            write(file, "dE", ΔE)
            write(file, "t_run", t_run)
            write(file, "mu_lims", μ_lims)
            write(file, "h_atm", z)
            write(file, "I0", I0)
            write(file, "mu_scatterings", scattering_data)
        close(file)
    end

    return nothing
end










# ======================================================================================== #
#                                  HEATING RATES                                         #
# ======================================================================================== #

"""
    make_heating_rate_file(directory_to_process)

Reads into a folder `directory_to_process` containing results from an AURORA.jl simulation,
loads the particle flux `Ie` (#e⁻/m²/s), and calculates the heating rate of thermal electrons
by superthermal electrons.

The heating rate is the rate at which energy is transferred from superthermal electrons to
thermal electrons through Coulomb collisions. It is saved as a function of altitude and time.

# Calling
`make_heating_rate_file(directory_to_process)`

# Inputs
- `directory_to_process`: absolute or relative path to the simulation directory to process.
"""
function make_heating_rate_file(directory_to_process)
    ## Find the files to process
    files = readdir(directory_to_process, join=true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]
    # The files are sorted in lexicographical order, so IeFlickering-100.mat will be loaded
    # before "IeFlickering-11.mat. We fix that with the following line which sorts them by
    # the number in the filename.
    sort!(files_to_process, by = x -> parse(Int, match(r"IeFlickering-(\d+)\.mat", basename(x))[1]))

    ## Load simulation grid
    f = matopen(files_to_process[1])
        z = read(f, "h_atm")
        E_centers = read(f, "E_centers")
    close(f)

    ## Load thermal electron density and temperature
    data = matread(joinpath(directory_to_process, "neutral_atm.mat"))
    ne = data["ne"]
    Te = data["Te"]

    ## Initialize arrays to store the results for each time-slice
    heating_rate = Vector{Matrix{Float64}}()
    t = Vector{}()

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
        n_z = length(z)
        n_μ = size(Ie_ztE, 1) ÷ n_z # (÷ returns an Int)
        n_t = size(Ie_ztE, 2)
        n_E = size(Ie_ztE, 3)
        Ie_ztE_omni = zeros(n_z, n_t, n_E)
        @views for i_μ in 1:n_μ
            Ie_ztE_omni .+= Ie_ztE[(i_μ - 1) * n_z .+ (1:n_z), :, :]
        end

        ## Calculate heating rate
        heating_rate_local = calculate_heating_rate(z, t_local, Ie_ztE_omni, E_centers, ne, Te)

        ## Push the heating rate of the current time-slice into a vector
        push!(heating_rate, heating_rate_local)
        push!(t, t_local)

        next!(p)
    end

    ## Concatenate along time
    heating_rate = reduce(hcat, heating_rate)
    t = reduce(vcat, t)

    ## Save results
    savefile = joinpath(directory_to_process, "heating_rate.mat")
    f = matopen(savefile, "w")
        write(f, "h_atm", z)
        write(f, "t", t)
        write(f, "heating_rate", heating_rate)
    close(f)

    println("Heating rates saved in $savefile")

    return nothing
end


"""
    calculate_heating_rate(z, t, Ie_ztE_omni, E_centers, ne, Te)

Calculate the heating rate of thermal electrons by superthermal electrons through Coulomb
collisions. The heating rate is the rate at which energy is transferred from superthermal
electrons to thermal electrons.

# Calling
`heating_rate = calculate_heating_rate(z, t, Ie_ztE_omni, E_centers, ne, Te)`

# Inputs
- `z`: altitude (m). Vector [n\\_z]
- `t`: time (s). Vector [n\\_t]
- `Ie_ztE_omni`: omnidirectional electron flux (#e⁻/m²/s). 3D array [n\\_z, n\\_t, n\\_E]
- `E_centers`: energy bin centers (eV). Vector [n\\_E]
- `ne`: thermal electron density (m⁻³). Vector [n\\_z]
- `Te`: thermal electron temperature (K). Vector [n\\_z]

# Output
- `heating_rate`: heating rate (eV/m³/s). 2D array [n\\_z, n\\_t]
"""
function calculate_heating_rate(z, t, Ie_ztE_omni, E_centers, ne, Te)
    # Initialize
    n_z = length(z)
    n_t = length(t)
    n_E = length(E_centers)
    heating_rate = zeros(n_z, n_t)

    # Calculate the energy loss rate to thermal electrons for each energy
    # loss_to_thermal_electrons returns the energy loss rate in eV/m
    L_th = zeros(n_z, n_E)
    for i_E in eachindex(E_centers)
        L_th[:, i_E] .= loss_to_thermal_electrons(E_centers[i_E], ne, Te)
    end

    # Calculate heating rate for each time step
    # The heating rate is the product of the flux and the energy loss rate, integrated over energy
    @views for i_t in eachindex(t)
        heating_rate[:, i_t] .= sum(Ie_ztE_omni[:, i_t, :] .* L_th, dims=2)
    end

    return heating_rate
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

"""
    make_Ie_top_file(sim::AuroraSimulation)

Convenience wrapper that calls [`make_Ie_top_file`](@ref) on `sim.savedir`.
"""
make_Ie_top_file(sim::AuroraSimulation) = make_Ie_top_file(sim.savedir)

"""
    make_current_file(sim::AuroraSimulation)

Convenience wrapper that calls [`make_current_file`](@ref) on `sim.savedir`.
"""
make_current_file(sim::AuroraSimulation) = make_current_file(sim.savedir)

"""
    make_heating_rate_file(sim::AuroraSimulation)

Convenience wrapper that calls [`make_heating_rate_file`](@ref) on `sim.savedir`.
"""
make_heating_rate_file(sim::AuroraSimulation) = make_heating_rate_file(sim.savedir)

"""
    downsampling_fluxes(sim::AuroraSimulation, downsampling_factor)

Convenience wrapper that calls [`downsampling_fluxes`](@ref) on `sim.savedir`.
"""
downsampling_fluxes(sim::AuroraSimulation, downsampling_factor) = downsampling_fluxes(sim.savedir, downsampling_factor)

"""
    make_psd_file(sim::AuroraSimulation; kwargs...)

Convenience wrapper that calls [`make_psd_file`](@ref) on `sim.savedir`.
"""
make_psd_file(sim::AuroraSimulation; kwargs...) = make_psd_file(sim.savedir; kwargs...)
