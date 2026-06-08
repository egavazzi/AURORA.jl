using NCDatasets: NCDataset, defDim, defVar


# ======================================================================================== #
#                              FLUX EXTRACTION & DOWNSAMPLING                            #
# ======================================================================================== #

"""
    make_Ie_top_file(directory_to_process)

Read `simulation_data.nc` from `directory_to_process`, extract the electron flux
at the top of the ionosphere (highest altitude), and write results to
`analysis/Ie_top.nc`.

# Calling
`make_Ie_top_file(directory_to_process)`

# Inputs
- `directory_to_process`: absolute or relative path to the simulation directory.
"""
function make_Ie_top_file(directory_to_process)
    ## Load simulation results
    result    = load_results(directory_to_process)
    Ie        = result.Ie          # [n_z, n_μ, n_t, n_E]
    t         = result.t
    E_centers = result.E_centers
    ΔE      = result.ΔE
    μ_lims  = result.μ_lims

    ## Beam weights
    θ_lims = acosd.(μ_lims)
    Ω_beam = beam_weight(θ_lims)   # [n_μ]

    ## Extract top altitude slice: Ie[end, :, :, :] → [n_μ, n_t, n_E]
    Ie_top_raw = Ie[end, :, :, :]      # in #e⁻/m²/s
    Ie_top = Ie_top_raw ./ reshape(ΔE, (1, 1, :)) ./ reshape(Ω_beam, (:, 1, 1))  # eV⁻¹ m⁻² s⁻¹ sr⁻¹

    ## Write to analysis/Ie_top.nc
    analysis_dir = joinpath(directory_to_process, "analysis")
    mkpath(analysis_dir)
    savefile = joinpath(analysis_dir, "Ie_top.nc")

    n_μ, n_t, n_E = size(Ie_top)
    NCDataset(savefile, "c") do ds
        defDim(ds, "pitch_angle",   n_μ)
        defDim(ds, "time",          n_t)
        defDim(ds, "energy",        n_E)
        defDim(ds, "energy_bounds",      n_E + 1)
        defDim(ds, "pitch_angle_bounds", n_μ + 1)

        pa_v = defVar(ds, "pitch_angle", Float64, ("pitch_angle",);
                      attrib=["units" => "1", "long_name" => "cosine of pitch angle"])
        pa_v[:] = mu_avg(θ_lims)

        t_v = defVar(ds, "time", Float64, ("time",);
                     attrib=["units" => "s", "long_name" => "simulation time"])
        t_v[:] = t

        en_v = defVar(ds, "energy", Float64, ("energy",);
                      attrib=["units" => "eV", "long_name" => "energy bin center"])
        en_v[:] = E_centers

        ee_v = defVar(ds, "energy_edges", Float64, ("energy_bounds",);
                      attrib=["units" => "eV", "long_name" => "energy bin edges"])
        ee_v[:] = vcat(E_centers[1] - ΔE[1]/2, E_centers .+ ΔE./2)

        ml_v = defVar(ds, "mu_lims", Float64, ("pitch_angle_bounds",);
                      attrib=["units" => "1", "long_name" => "pitch-angle cosine bin boundaries"])
        ml_v[:] = μ_lims

        bw_v = defVar(ds, "beam_weight", Float64, ("pitch_angle",);
                      attrib=["units" => "sr", "long_name" => "solid-angle beam weight"])
        bw_v[:] = Ω_beam

        Ie_top_raw_v = defVar(ds, "Ie_top_raw", Float32, ("pitch_angle", "time", "energy");
                              deflatelevel=4,
                              attrib=["units"     => "m-2 s-1",
                                       "long_name" => "electron flux at top of ionosphere"])
        Ie_top_raw_v[:, :, :] = Float32.(Ie_top_raw)

        Ie_top_v = defVar(ds, "Ie_top", Float32, ("pitch_angle", "time", "energy");
                          deflatelevel=4,
                          attrib=["units"     => "eV-1 m-2 s-1 sr-1",
                                   "long_name" => "differential electron flux at top of ionosphere"])
        Ie_top_v[:, :, :] = Float32.(Ie_top)
    end

    println("Top flux saved in $savefile")
    return nothing
end


# ======================================================================================== #
#                                  FIELD-ALIGNED CURRENTS                                #
# ======================================================================================== #

"""
    make_current_file(directory_to_process)

Read `simulation_data.nc` from `directory_to_process`, compute field-aligned
current density and energy flux, and write results to `analysis/currents.nc`.

# Calling
`make_current_file(directory_to_process)`
"""
function make_current_file(directory_to_process)
    ## Load simulation results
    result    = load_results(directory_to_process)
    Ie        = result.Ie          # [n_z, n_μ, n_t, n_E]
    t         = result.t
    z         = result.h_atm
    E_centers = result.E_centers
    μ_lims    = result.μ_lims

    θ_lims   = acosd.(μ_lims)
    μ_center = mu_avg(θ_lims)

    ## Elementary charge
    q_e = 1.602176620898e-19  # C

    n_z, n_μ, n_t, n_E = size(Ie)
    J_up      = zeros(n_z, n_t)
    J_down    = zeros(n_z, n_t)
    IeE_up    = zeros(n_z, n_t)
    IeE_down  = zeros(n_z, n_t)

    @views for i_μ in 1:n_μ
        if μ_center[i_μ] > 0
            J_up    .+= q_e * abs(μ_center[i_μ]) .* dropdims(sum(Ie[:, i_μ, :, :], dims=3), dims=3)
            IeE_up  .+= abs(μ_center[i_μ]) .* dropdims(sum(Ie[:, i_μ, :, :] .* reshape(E_centers, 1, 1, :), dims=3), dims=3)
        else
            J_down   .+= q_e * abs(μ_center[i_μ]) .* dropdims(sum(Ie[:, i_μ, :, :], dims=3), dims=3)
            IeE_down .+= abs(μ_center[i_μ]) .* dropdims(sum(Ie[:, i_μ, :, :] .* reshape(E_centers, 1, 1, :), dims=3), dims=3)
        end
    end

    ## Write to analysis/currents.nc
    analysis_dir = joinpath(directory_to_process, "analysis")
    mkpath(analysis_dir)
    savefile = joinpath(analysis_dir, "currents.nc")

    NCDataset(savefile, "c") do ds
        n_z = length(z)
        n_t = length(t)
        defDim(ds, "altitude", n_z)
        defDim(ds, "time",     n_t)

        alt_v = defVar(ds, "altitude", Float64, ("altitude",);
                       attrib=["units" => "m", "long_name" => "altitude"])
        alt_v[:] = z
        t_v = defVar(ds, "time", Float64, ("time",);
                     attrib=["units" => "s", "long_name" => "simulation time"])
        t_v[:] = t

        for (name, data, long_name, units) in (
                ("J_up",     J_up,     "Upward field-aligned current density",  "A m-2"),
                ("J_down",   J_down,   "Downward field-aligned current density", "A m-2"),
                ("IeE_up",   IeE_up,   "Upward field-aligned energy flux",       "eV m-2 s-1"),
                ("IeE_down", IeE_down, "Downward field-aligned energy flux",     "eV m-2 s-1"),
            )
            v = defVar(ds, name, Float64, ("altitude", "time");
                       deflatelevel=4,
                       attrib=["units" => units, "long_name" => long_name])
            v[:, :] = data
        end
    end

    println("Currents saved in $savefile")
    return nothing
end


# ======================================================================================== #
#                         AuroraSimulation convenience wrappers                          #
# ======================================================================================== #

"""
    make_Ie_top_file(sim::AuroraSimulation)

Convenience wrapper that calls [`make_Ie_top_file`](@ref) on `sim.output.savedir`.
"""
make_Ie_top_file(sim::AuroraSimulation) = make_Ie_top_file(sim.output.savedir)

"""
    make_current_file(sim::AuroraSimulation)

Convenience wrapper that calls [`make_current_file`](@ref) on `sim.output.savedir`.
"""
make_current_file(sim::AuroraSimulation) = make_current_file(sim.output.savedir)
