using MAT
using Interpolations

# conservation of fzvzmu and Ie ✅
function convert_M_to_I(path_to_vlasov_initial_file, path_to_vlasov_simulation,
    index_specie, HMR_VZ, HMR_MU, E_max, θ_lims, first_run=0)
    ## Import the data
    B, dB, Nz, z, zcorn, dz, dt, particle = load_Bfield(path_to_vlasov_simulation);

    dt, Niter, dump_period_fields, fields_per_file, dump_period_distr, dump_period_distr_IonoBoundary,
    dump_period_distr_1v, dump_start, shift_test_period, resistance, Nz, zmin, zmax,
    Nspecies, const_a, BC_Poisson, voltage, initialiser, voltage_init, E0, startfromdumpfile,
    transffilename, voltagefilename, _ = read_input(joinpath(path_to_vlasov_simulation, "inputb6.m"));

    # Extract Nvz, Nmu and m
    Nvz = Int64(particle[index_specie].Nvz)
    Nmu = Int64(particle[index_specie].Nmu)
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
    Bz = B[index_z]
    for i_t in 1:Nt # Loop over time
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

    # Calculate the coordinates in the (E, μ_pitch)-grid of all the points in the (vz, μ_mag)-finer grid
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




## -------------------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------------------- ##
## -------------------------------------------------------------------------------------- ##




using ProgressMeter
using QuadGK

# speedup ~240x
# conservation of Ie and fzvzmu ✅
function convert_I_to_M(path_to_vlasov_simulation, path_to_aurora_simulation,
    index_specie, HMR_MU, HMR_E, E_max, θ_lims)
    ## Import the data
    B, dB, Nz, z, zcorn, dz, dt, particle = load_Bfield(path_to_vlasov_simulation);

    dt, Niter, dump_period_fields, fields_per_file, dump_period_distr, dump_period_distr_IonoBoundary,
    dump_period_distr_1v, dump_start, shift_test_period, resistance, Nz, zmin, zmax,
    Nspecies, const_a, BC_Poisson, voltage, initialiser, voltage_init, E0, startfromdumpfile,
    transffilename, voltagefilename, _ = read_input(joinpath(path_to_vlasov_simulation, "inputb6.m"));

    # Extract Nvz, Nmu, m and B
    Nvz = Int64(particle[index_specie].Nvz)
    Nmu = Int64(particle[index_specie].Nmu)
    m = particle[index_specie].mass
    Bz = B[Nz - 1] # field at one point over the ionosphere (but remember that they overlap)
    ## Extracting Ie
    file = matopen(joinpath(path_to_aurora_simulation, "Ie_MI_top.mat"))
        Ie1 = read(file, "Ie_top_raw")
        Ie2 = read(file, "Ie_top2_raw")
        E_grid = vec(read(file, "E"))
    close(file)
    ## Extract the (vz, μ_mag)-grid
    vz = particle[index_specie].vz
    Δvz = particle[index_specie].dvz
    μ_mag = particle[index_specie].mu
    Δμ_mag = particle[index_specie].dmu
    ## Resize E-grid 
    iE_max = findmin(abs.(E_grid .- E_max))[2];  # find the index for the upper limit of the energy grid
    E_grid = E_grid[1:iE_max]                   # crop E accordingly
    ΔE = diff(E_grid); ΔE = [ΔE; ΔE[end]]
    ## Refine E-grid (using gridded interpolation)
    F = interpolate((eachindex(E_grid), ), E_grid, Gridded(Linear()))
    F = extrapolate(F, Line()) # to be able to take one point over the max of E
    E_grid_finer = F(1:(1/HMR_E):length(E_grid)+1) # ⬅ there
    ΔE_finer = 1 / HMR_E .* ΔE
    ΔE_finer = repeat(ΔE_finer, inner=HMR_E)
    E_finer = E_grid_finer[2:end] - 0.5 * ΔE_finer
    ## Refine μ_pitch-grid (using gridded interpolation)
    μ_pitch_grid = cosd.(θ_lims)
    G = interpolate((eachindex(θ_lims), ), θ_lims, Gridded(Linear()))
    θ_lims_finer = G(1:(1/HMR_MU):length(θ_lims))
    μ_pitch_grid_finer = cosd.(θ_lims_finer)
    μ_pitch_finer = mu_avg(θ_lims_finer)
    BeamW = Vector{Float64}(undef, length(μ_pitch_grid) - 1)
    for i_μ in eachindex(BeamW) #(length(μ_pitch_grid)-1):-1:1
        BeamW[i_μ] = 2 * π * abs(quadgk(sin, acos(μ_pitch_grid[i_μ]), acos(μ_pitch_grid[i_μ + 1]))[1])
    end
    BeamW_finer = Vector{Float64}(undef, length(μ_pitch_grid_finer) - 1)
    for i_μ in eachindex(BeamW_finer) #(length(μ_pitch_grid_finer)-1):-1:1
        BeamW_finer[i_μ] = 2 * π * abs(quadgk(sin, acos(μ_pitch_grid_finer[i_μ]), acos(μ_pitch_grid_finer[i_μ + 1]))[1])
    end

    Nt = size(Ie1, 2)
    f_out1 = zeros(Nt, Nvz, Nmu)
    f_out2 = zeros(Nt, Nvz, Nmu)
    p = Progress(Nt, desc=string("Converting the e⁻ fluxes"))
    Threads.@threads for i_t in 1:Nt # Loop over time
        f_out1[i_t, :, :] = convert_Ie_to_fzvzmu(E_finer, μ_pitch_finer, BeamW, BeamW_finer, vz, μ_mag, Δvz, Δμ_mag, 
                                        Bz, m, HMR_MU, HMR_E, Ie1[:, i_t, :])
        f_out2[i_t, :, :] = convert_Ie_to_fzvzmu(E_finer, μ_pitch_finer, BeamW, BeamW_finer, vz, μ_mag, Δvz, Δμ_mag, 
                                        Bz, m, HMR_MU, HMR_E, Ie2[:, i_t, :])
        next!(p)
    end

    # Reshape f from [Nt, Nvz, Nmu] to [Nt * Nvz, Nmu]
    # ➡ we stack the [Nvz, Nmu] matrices vertically
    f_out1 = reshape(permutedims(f_out1, (2, 1, 3)), (:, size(f_out1, 3)));
    f_out2 = reshape(permutedims(f_out2, (2, 1, 3)), (:, size(f_out2, 3)));

    return f_out1, f_out2
end


"""
    convert_Ie_to_fzvzmu(E_finer, μ_pitch_finer, BeamWeight, BeamWeight_finer, vz, μ_mag, Δvz, Δμ_mag, 
    Bz, m, HMR_MU, HMR_E, Ie)

This function converts a particle flux Ie (#e⁻/m²/s) into a function distribution fzvzmu (#e⁻/m⁶/s³)

The function distribution fzvzmu is defined along a magnetic field line and over a (v\\_z, magnetic\\_moment)-grid
The particle flux Ie is defined along a magnetic field line and over an (Energy, pitch_angle)-grid

# Calling
`fzvzmu = convert_Ie_to_fzvzmu(E_finer, μ_pitch_finer, BeamWeight, BeamWeight_finer, vz, μ_mag, Δvz, Δμ_mag, 
Bz, m, HMR_MU, HMR_E, Ie)`
"""
function convert_Ie_to_fzvzmu(E_finer, μ_pitch_finer, BeamWeight, BeamWeight_finer, vz, μ_mag, Δvz, Δμ_mag, 
    Bz, m, HMR_MU, HMR_E, Ie)
    
    fzvzmu = zeros(length(vz), length(μ_mag))
    Ie_finer = Ie ./ (BeamWeight ./ sum(BeamWeight))
    Ie_finer = repeat(Ie_finer ./ HMR_E, inner=(HMR_MU, HMR_E))
    Ie_finer = Ie_finer .* (BeamWeight_finer ./ sum(BeamWeight_finer))

    # Calculate the coordinate in the (vz, μ_mag)-grid of all the points in the (E, μ_pitch)-finer grid
    vz_local = [- sqrt(2 * Y * 1.6e-19 / m) * X for X in μ_pitch_finer, Y in E_finer]
    μ_mag_local = [Y * 1.6e-19 / Bz * (1 - X^2) for X in μ_pitch_finer, Y in E_finer]
    idx_vz = findnearestindex.(Ref(vz), vz_local)
    idx_μ_mag = findnearestindex.(Ref(μ_mag), μ_mag_local)

    v = [sqrt(i^2 + 2 * Bz / m * j) for i in vz, j in μ_mag]

    for i in eachindex(μ_pitch_finer)
        for j in eachindex(E_finer)
            # convert Ie into fzvzmu
            fzvzmu[idx_vz[i, j], idx_μ_mag[i, j]] += 1 / v[idx_vz[i, j], idx_μ_mag[i, j]] *
                                                        Ie_finer[i, j] / Δvz / Δμ_mag[idx_μ_mag[i, j]]
        end
    end
    return fzvzmu
end


