using NCDatasets: NCDataset, defDim, defVar


# ======================================================================================== #
#                                  HEATING RATES                                         #
# ======================================================================================== #

"""
    make_heating_rate_file(directory_to_process)

Read `simulation_data.nc` and `inputs/atmosphere.nc` from `directory_to_process`,
compute the thermal-electron heating rate by superthermal electrons, and write
results to `analysis/heating_rate.nc`.

# Calling
`make_heating_rate_file(directory_to_process)`

# Inputs
- `directory_to_process`: absolute or relative path to the simulation directory.
"""
function make_heating_rate_file(directory_to_process)
    ## Load simulation results
    result    = load_results(directory_to_process)
    Ie        = result.Ie          # [n_z, n_μ, n_t, n_E]
    t         = result.t
    z         = result.h_atm
    E_centers = result.E_centers

    ## Load thermal electron density and temperature
    atm = read_atmosphere_nc(directory_to_process)
    ne  = atm.ne
    Te  = atm.Te

    ## Sum Ie over pitch-angle beams → omnidirectional flux [n_z, n_t, n_E]
    Ie_omni = dropdims(sum(Ie, dims=2), dims=2)

    ## Calculate heating rate
    heating_rate = calculate_heating_rate(z, t, Ie_omni, E_centers, ne, Te)

    ## Write to analysis/heating_rate.nc
    analysis_dir = joinpath(directory_to_process, "analysis")
    mkpath(analysis_dir)
    savefile = joinpath(analysis_dir, "heating_rate.nc")

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

        hr_v = defVar(ds, "heating_rate", Float64, ("altitude", "time");
                      deflatelevel=4,
                      attrib=["units"     => "eV m-3 s-1",
                               "long_name" => "thermal electron heating rate by superthermal electrons"])
        hr_v[:, :] = heating_rate
    end

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
    n_z = length(z)
    n_t = length(t)
    n_E = length(E_centers)
    heating_rate = zeros(n_z, n_t)

    # Energy loss rate to thermal electrons per unit path length (eV/m), for each energy bin
    L_th = zeros(n_z, n_E)
    for i_E in eachindex(E_centers)
        L_th[:, i_E] .= loss_to_thermal_electrons(E_centers[i_E], ne, Te)
    end

    # Heating rate = flux x energy-loss-rate, integrated over energy
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

Convenience wrapper that calls [`make_heating_rate_file`](@ref) on `sim.output.savedir`.
"""
make_heating_rate_file(sim::AuroraSimulation) = make_heating_rate_file(sim.output.savedir)
