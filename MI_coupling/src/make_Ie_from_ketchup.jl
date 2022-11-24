using MAT
using Interpolations
using Aurora
using ProgressMeter
"""
    make_Ie_from_ketchup(path_to_vlasov_initial_file, path_to_vlasov_simulation,
        index_specie, zmax_to_extract, h_atm_iono, HMR_VZ, HMR_MU,
        iE_max, θ_lims)

This function extracts and converts the function distributions from ketchup into fluxes of
electrons Ie.

It looks through a directory with all the fzvzmu***.mat files from a simulation. Each file 
should correspond to one time step and contain the function distribution over all the space 
points of the simulation.

Speedup compared to the Matlab version : ~100x

# Calling
`Ie, z_DL =  make_Ie_from_ketchup(path_to_vlasov_initial_file, path_to_vlasov_simulation,
index_specie, zmax_to_extract, h_atm_iono, HMR_VZ, HMR_MU,
iE_max, θ_lims)`
"""
function make_Ie_from_ketchup(path_to_vlasov_initial_file, path_to_vlasov_simulation,
                                index_specie, zmax_to_extract, h_atm_iono, HMR_VZ, HMR_MU,
                                E_max, θ_lims)

    ## Import the data
    B, dB, Nz, z, zcorn, dz, dt, particle = load_Bfield(path_to_vlasov_simulation);

    dt, Niter, dump_period_fields, fields_per_file, dump_period_distr, dump_period_distr_IonoBoundary,
    dump_period_distr_1v, dump_start, shift_test_period, resistance, Nz, zmin, zmax,
    Nspecies, const_a, BC_Poisson, voltage, initialiser, voltage_init, E0, startfromdumpfile,
    transffilename, voltagefilename, _ = read_input(joinpath(path_to_vlasov_simulation, "inputb6.m"));

    # Extract Nvz, Nmu and m
    Nvz = particle[index_specie].Nvz
    Nmu = particle[index_specie].Nmu
    m = particle[index_specie].mass
    ## Extracting fzvzmu_in
    # fzvzmu = load_fzvzmu_parallel(path_to_vlasov_initial_file, path_to_vlasov_simulation, index_specie, 6);
    fzvzmu = load_fzvzmu_serial(path_to_vlasov_initial_file, path_to_vlasov_simulation, index_specie);
    ## Create z_DL : altitude (in km) above the ionosphere
    z_DL = z[end] .- z[end:-1:1] .+ h_atm_iono[end - 1]
    z_DL = z_DL ./ 1e3
    i_zmax = findmin(abs.(z_DL .- zmax_to_extract))[2]
    z_DL = z_DL[1:i_zmax]
    ## Extract the (vz, μ_mag)-grid
    vz = particle[index_specie].vz
    vz_grid = particle[index_specie].vzcorn
    Δvz = particle[index_specie].dvz
    μ_mag = particle[index_specie].mu
    μ_mag_grid = particle[index_specie].mucorn
    Δμ_mag = particle[index_specie].dmu
    ## Refine vz_grid (using gridded interpolation)
    F = interpolate((eachindex(vz_grid), ), vz_grid, Gridded(Linear()))
    vz_grid_finer = F(1:1/HMR_VZ:length(vz_grid))
    Δvz_finer = Δvz * 1/HMR_VZ
    vz_finer = vz_grid_finer[2:end] .- 0.5 * Δvz_finer
    ## Refine μ_mag_grid (using gridded interpolation)
    G = interpolate((eachindex(μ_mag_grid), ), μ_mag_grid, Gridded(Linear()))
    μ_mag_grid_finer = G(1:1/HMR_MU:length(μ_mag_grid))
    Δμ_mag_finer = 1/HMR_MU .* Δμ_mag; 
    Δμ_mag_finer = repeat(Δμ_mag_finer, inner=HMR_MU)
    μ_mag_finer = μ_mag_grid_finer[2:end] - 0.5 * Δμ_mag_finer

    ## Make E-grid
    E_function(X, dE_initial, dE_final, C, X0) = dE_initial + (1 + tanh(C * (X - X0))) / 2 * dE_final
    E = cumsum(E_function.(0:2000, 0.15, 11.5, 0.05, 80)) .+ 1.9
    iE_max = findmin(abs.(E .- E_max))[2];  # find the index for the upper limit of the energy grid
    E = E[1:iE_max];                        # crop E accordingly
    dE = diff(E); dE = [dE; dE[end]]
    ## Make μ_pitch-grid
    μ_pitch = mu_avg(θ_lims)

    Nt = size(fzvzmu, 1)
    Ie = zeros(i_zmax * length(μ_pitch), Nt, iE_max)
    p = Progress(i_zmax, desc=string("Converting the function distribution"))
    println("Number of altitudes : ", i_zmax)
    println("Number of timesteps : ", Nt)
    Threads.@threads for i_z in 1:i_zmax # loop over altitude
        i_zz = Nz - i_z + 1 # because ketchup matrices are flipped
        Bz = B[i_zz]
        for i_t in 1:Nt # loop over time
            Ie_temp = convert_fzvzmu_to_Ie(vz, μ_mag, vz_finer, μ_mag_finer, Δvz_finer, Δμ_mag_finer, 
                                            E, μ_pitch, Bz, m, HMR_VZ, HMR_MU, fzvzmu[i_t, :, :, i_zz])
            for i_μ in 1:length(μ_pitch)
                Ie[i_zmax * (i_μ - 1) + i_z, i_t, :] = Ie_temp[:, i_μ]
            end
        end
        next!(p)
    end
    # something to save the data?
    return Ie, z_DL
end




function convert_M_to_I(path_to_vlasov_initial_file, path_to_vlasov_simulation,
                        index_specie, HMR_VZ, HMR_MU, E_max, θ_lims, first_run=0)
    ## Import the data
    # should probably clean this comment block ⬇⬇⬇
    # if first_run != 0
    #     # then in that case there is no Bfield.mat in the vlasov_simulation folder as it has
    #     # not been generated yet, so we need to load B from the initial file folder instead.
    #     idx_folder = findall("/", path_to_vlasov_initial_file)[end - 1][1]
    #     path_to_initial_folder = path_to_vlasov_initial_file[1:idx_folder]
    #     B, dB, Nz, z, zcorn, dz, dt, particle = load_Bfield(path_to_initial_folder);
    # else
    B, dB, Nz, z, zcorn, dz, dt, particle = load_Bfield(path_to_vlasov_simulation);
    # end

    dt, Niter, dump_period_fields, fields_per_file, dump_period_distr, dump_period_distr_IonoBoundary,
    dump_period_distr_1v, dump_start, shift_test_period, resistance, Nz, zmin, zmax,
    Nspecies, const_a, BC_Poisson, voltage, initialiser, voltage_init, E0, startfromdumpfile,
    transffilename, voltagefilename, _ = read_input(joinpath(path_to_vlasov_simulation, "inputb6.m"));

    # Extract Nvz, Nmu and m
    Nvz = particle[index_specie].Nvz
    Nmu = particle[index_specie].Nmu
    m = particle[index_specie].mass
    ## Extracting fzvzmu_in
    fzvzmu = load_fzvzmuIB_serial(path_to_vlasov_initial_file, path_to_vlasov_simulation, index_specie, first_run);
    ## Extract the (vz, μ_mag)-grid
    vz = particle[index_specie].vz
    vz_grid = particle[index_specie].vzcorn
    Δvz = particle[index_specie].dvz
    μ_mag = particle[index_specie].mu
    μ_mag_grid = particle[index_specie].mucorn
    Δμ_mag = particle[index_specie].dmu
    ## Refine vz_grid (using gridded interpolation)
    F = interpolate((eachindex(vz_grid), ), vz_grid, Gridded(Linear()))
    vz_grid_finer = F(1:1/HMR_VZ:length(vz_grid))
    Δvz_finer = Δvz * 1/HMR_VZ
    vz_finer = vz_grid_finer[2:end] .- 0.5 * Δvz_finer
    ## Refine μ_mag_grid (using gridded interpolation)
    G = interpolate((eachindex(μ_mag_grid), ), μ_mag_grid, Gridded(Linear()))
    μ_mag_grid_finer = G(1:1/HMR_MU:length(μ_mag_grid))
    Δμ_mag_finer = 1/HMR_MU .* Δμ_mag; 
    Δμ_mag_finer = repeat(Δμ_mag_finer, inner=HMR_MU)
    μ_mag_finer = μ_mag_grid_finer[2:end] - 0.5 * Δμ_mag_finer
    ## Make E-grid
    E_function(X, dE_initial, dE_final, C, X0) = dE_initial + (1 + tanh(C * (X - X0))) / 2 * dE_final
    E = cumsum(E_function.(0:2000, 0.15, 11.5, 0.05, 80)) .+ 1.9
    iE_max = findmin(abs.(E .- E_max))[2];  # find the index for the upper limit of the energy grid
    E = E[1:iE_max];                        # crop E accordingly
    dE = diff(E); dE = [dE; dE[end]]
    ## Make μ_pitch-grid
    μ_pitch = mu_avg(θ_lims)

    index_z = Nz - 1; # we extract the function distribution at one point over the ionospheric boundary
    Nt = size(fzvzmu, 1)
    Ie = zeros(length(μ_pitch), Nt, iE_max)
    #  p = Progress(Nt, desc=string("Converting the function distribution"))
    #  println("Number of timesteps : ", Nt)
    Bz = B[index_z]
    for i_t in 1:Nt # loop over time
        Ie_temp = convert_fzvzmu_to_Ie(vz, μ_mag, vz_finer, μ_mag_finer, Δvz_finer, Δμ_mag_finer, 
                                        E, μ_pitch, Bz, m, HMR_VZ, HMR_MU, fzvzmu[i_t, :, :, 1])
        for i_μ in 1:length(μ_pitch)
            Ie[i_μ, i_t, :] = Ie_temp[:, i_μ]
        end
    end
    #  next!(p)
    return Ie
end




"""
    function convert_fzvzmu_to_Ie(vz, μ_mag, vz_finer, μ_mag_finer, Δvz_finer, Δμ_mag_finer, 
        E, μ_pitch, Bz, m, HMR_VZ, HMR_MU, fzvzmu)

This function converts a particle function distribution fzvzmu (#e⁻/m⁶/s³) into a particle flux Ie (#e⁻/m²/s)

The function distribution fzvzmu is defined along a magnetic field line and over a (v\\_z, magnetic\\_moment)-grid
The particle flux Ie is defined along a magnetic field line and over an (Energy, pitch_angle)-grid

# Calling
`Ie = convert_fzvzmu_to_Ie(vz, μ_mag, vz_finer, μ_mag_finer, Δvz_finer, Δμ_mag_finer, 
E, μ_pitch, Bz, m, HMR_VZ, HMR_MU, fzvzmu)`
"""
function convert_fzvzmu_to_Ie(vz, μ_mag, vz_finer, μ_mag_finer, Δvz_finer, Δμ_mag_finer, 
                                E, μ_pitch, Bz, m, HMR_VZ, HMR_MU, fzvzmu)
    Ie = zeros(length(E), length(μ_pitch))
    fzvzmu_finer = repeat(fzvzmu, inner=(HMR_VZ, HMR_MU))
    
    v = [sqrt(i^2 + 2 * Bz / m * j) for i in vz, j in μ_mag]
    v = repeat(v, inner=(HMR_VZ, HMR_MU))

    # Calculate the coordinates in the (E, μ_pitch)-grid for all the points in the (vz, μ_mag)-grid
    μ_pitch_local = [-sign(X) * cos(atan(sqrt(2 * Bz / m * Y) / X)) for X in vz_finer, Y in μ_mag_finer]
    E_local = [0.5 * m * (X^2 + 2 * Bz / m * Y) / 1.6e-19 for X in vz_finer, Y in μ_mag_finer]
    idx_pitch = findnearestindex.(Ref(μ_pitch), μ_pitch_local) 
    idx_energy = findnearestindex.(Ref(E), E_local)
    
    for i in eachindex(vz_finer)
        for j in eachindex(μ_mag_finer)
            # convert fzvzmu into Ie
            Ie[idx_energy[i, j], idx_pitch[i, j]] += v[i, j] * fzvzmu_finer[i, j] * Δvz_finer * Δμ_mag_finer[j]
        end
    end
    return Ie
end