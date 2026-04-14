using MAT: matopen, matread
using ProgressMeter: Progress, next!


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

    if isempty(files_to_process)
        @warn "No simulation results found in $directory_to_process. Skipping heating rate calculations."
        return nothing
    end

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
        append!(t, t_local)

        next!(p)
    end

    ## Concatenate along time
    heating_rate = reduce(hcat, heating_rate)

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
        heating_rate[:, i_t] .= dropdims(sum(Ie_ztE_omni[:, i_t, :] .* L_th, dims=2); dims=2)
    end

    return heating_rate
end


# ======================================================================================== #
#                         AuroraSimulation convenience wrapper                           #
# ======================================================================================== #

"""
    make_heating_rate_file(sim::AuroraSimulation)

Convenience wrapper that calls [`make_heating_rate_file`](@ref) on `sim.savedir`.
"""
make_heating_rate_file(sim::AuroraSimulation) = make_heating_rate_file(sim.savedir)
