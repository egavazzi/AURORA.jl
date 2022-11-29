using MAT
using Interpolations
using QuadGK
using ProgressMeter


# conservation of Ie and fzvzvperp ✅
function make_f_from_Aurora(path_to_file, HMR_E, HMR_MU)
    m_e = 9.109e-31 # electron mass (kg)
    ## Extracting Ie
    file = matopen(path_to_file)
        Ie = read(file, "Ie_ztE") # [n_μ * nz, nt, nE]
        E = read(file, "E")
        μ_pitch_grid = vec(read(file, "mu_lims"))
        t_run = read(file, "t_run")
        h_atm = read(file, "h_atm")
    close(file)
    ## Resize E-grid
    E_grid = E[1:size(Ie, 3)]
    ΔE = diff(E_grid); ΔE = [ΔE; ΔE[end]]
    ## Refine E-grid (using gridded interpolation)
    F = interpolate((eachindex(E_grid), ), E_grid, Gridded(Linear()))
    F = extrapolate(F, Line()) # to be able to take one point over the max of E
    E_grid_finer = F(1:(1/HMR_E):length(E_grid)+1) # ⬅ there
    ΔE_finer = 1 / HMR_E .* ΔE
    ΔE_finer = repeat(ΔE_finer, inner=HMR_E)
    E_finer = E_grid_finer[2:end] .- 0.5 * ΔE_finer
    ## Refine μ_pitch-grid (using gridded interpolation)
    θ_lims = acosd.(μ_pitch_grid)
    G = interpolate((eachindex(θ_lims), ), θ_lims, Gridded(Linear()))
    θ_lims_finer = G(1:(1/HMR_MU):length(θ_lims))
    μ_pitch_grid_finer = cosd.(θ_lims_finer)
    μ_pitch_finer = mu_avg(θ_lims_finer)
    BeamW = Vector{Float64}(undef, length(μ_pitch_grid) - 1)
    for i_μ in eachindex(BeamW)
        BeamW[i_μ] = 2 * π * abs(quadgk(sin, acos(μ_pitch_grid[i_μ]), acos(μ_pitch_grid[i_μ + 1]))[1])
    end
    BeamW_finer = Vector{Float64}(undef, length(μ_pitch_grid_finer) - 1)
    for i_μ in eachindex(BeamW_finer)
        BeamW_finer[i_μ] = 2 * π * abs(quadgk(sin, acos(μ_pitch_grid_finer[i_μ]), acos(μ_pitch_grid_finer[i_μ + 1]))[1])
    end
    ## Make the vz grid
    vz_grid_plus = [sqrt(2 * X * 1.6e-19 / m_e) for X in E_grid]
    Δvz_plus = diff(vz_grid_plus); Δvz_plus = [Δvz_plus; Δvz_plus[end]]
    vz_plus = vz_grid_plus .+ Δvz_plus/2
    vz_minus = - reverse(vz_plus)
    Δvz_minus = reverse(Δvz_plus) # Δvz is always > 0
    vz = vcat(vz_minus, vz_plus) # vz takes negative and positive values
    Δvz = vcat(Δvz_minus, Δvz_plus)
    ## Make the v⟂ grid
    vperp_grid = [sqrt(2 * X * 1.6e-19 / m_e) for X in E_grid]
    Δvperp = diff(vperp_grid); Δvperp = [Δvperp; Δvperp[end]]
    vperp = vperp_grid .+ Δvz_minus/2 # whereas v⟂ takes only positive values

    Nz = length(h_atm)
    Nt = length(t_run)
    fzvzvperp = zeros(Nt, Nz, length(vz), length(vperp))
    p = Progress(Nz, desc=string("Converting Ie"))
    Threads.@threads for i_z in 1:Nz
        idx_z = [i_z + (i_μ - 1) * Nz for i_μ in 1:(length(θ_lims) - 1)]
        for i_t in 1:Nt
            fzvzvperp[i_t, i_z, :, :] = convert_Ie_to_fzvzvperp(E_finer, μ_pitch_finer,
                BeamW, BeamW_finer, vz, vperp, Δvz, Δvperp, m_e, HMR_MU, HMR_E, Ie[idx_z, i_t, :])
        end
        next!(p)
    end
    return fzvzvperp, E_grid, t_run, h_atm, μ_pitch_grid , vz, Δvz, vperp, Δvperp
end




function convert_Ie_to_fzvzvperp(E_finer, μ_pitch_finer, BeamWeight, BeamWeight_finer, vz, vperp, Δvz, Δvperp,
    m, HMR_MU, HMR_E, Ie)

    fzvzvperp = zeros(length(vz), length(vperp))
    Ie_finer = Ie ./ (BeamWeight ./ sum(BeamWeight))
    Ie_finer = repeat(Ie_finer ./ HMR_E, inner=(HMR_MU, HMR_E))
    Ie_finer = Ie_finer .* (BeamWeight_finer ./ sum(BeamWeight_finer))

    # Calculate the coordinate in the (vz, v⟂)-grid of all the points in the (E, μ_pitch)-finer grid
    vz_local = [- sqrt(2 * Y * 1.6e-19 / m) * X for X in μ_pitch_finer, Y in E_finer]
    vperp_local = [sqrt(2 * Y * 1.6e-19 / m) * sqrt(1 - X^2) for X in μ_pitch_finer, Y in E_finer]
    idx_vz = findnearestindex.(Ref(vz), vz_local)
    idx_vperp = findnearestindex.(Ref(vperp), vperp_local)

    v = [sqrt(i^2 + j^2) for i in vz, j in vperp]

    for i in eachindex(μ_pitch_finer)
        for j in eachindex(E_finer)
            # convert Ie into fzvzperp
            fzvzvperp[idx_vz[i, j], idx_vperp[i, j]] += 1 / v[idx_vz[i, j], idx_vperp[i, j]] *
                                                        Ie_finer[i, j] /
                                                        (Δvz[idx_vz[i, j]] * Δvperp[idx_vperp[i, j]]^2)
        end
    end
    return fzvzvperp
end
