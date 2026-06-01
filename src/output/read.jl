using NCDatasets: NCDataset

# ======================================================================================== #
#                              RESULT FILE READERS                                       #
# ======================================================================================== #

"""
    read_simulation_nc(sim_dir) → NamedTuple

Read `simulation_data.nc` from `sim_dir` and return a named tuple with all
simulation output arrays.

# Returns
Named tuple with fields:
- `Ie`        : electron number flux `[n_z, n_μ, n_t, n_E]` (m⁻² s⁻¹)
- `t`         : time axis `[n_t]` (seconds)
- `h_atm`     : altitude grid `[n_z]` (metres)
- `E_centers` : energy bin centres `[n_E]` (eV)
- `E_edges`   : energy bin edges `[n_E+1]` (eV)
- `dE`        : energy bin widths `[n_E]` (eV) (= diff of E_edges)
- `mu_lims`   : pitch-angle cosine bin boundaries `[n_μ+1]`

# Example
```julia
result = read_simulation_nc("my_run/")
result.t       # time axis (s)
result.Ie      # flux [n_z, n_μ, n_t, n_E]
```
"""
function read_simulation_nc(sim_dir::AbstractString)
    nc_path = joinpath(sim_dir, "simulation_data.nc")
    NCDataset(nc_path, "r") do ds
        Ie        = Float64.(Array(ds["Ie"]))    # [n_z, n_μ, n_t, n_E]
        t         = Array(ds["time"])
        h_atm     = Array(ds["altitude"])
        E_centers = Array(ds["energy"])
        E_edges   = Array(ds["energy_edges"])
        dE        = diff(E_edges)
        mu_lims   = Array(ds["mu_lims"])
        return (; Ie, t, h_atm, E_centers, E_edges, dE, mu_lims)
    end
end


"""
    read_atmosphere_nc(sim_dir) → NamedTuple

Read `inputs/atmosphere.nc` from `sim_dir`. Returns a named tuple with fields
`h_atm`, `ne`, `Te`, and one field per species (e.g. `nN2`, `nO2`, `nO`).
"""
function read_atmosphere_nc(sim_dir::AbstractString)
    nc_path = joinpath(sim_dir, "inputs", "atmosphere.nc")
    NCDataset(nc_path, "r") do ds
        h_atm = Array(ds["altitude"])
        ne    = Array(ds["ne"])
        Te    = Array(ds["Te"])
        # Collect all species density variables (variables starting with 'n' other than 'ne')
        species_vars = filter(k -> startswith(k, "n") && k != "ne" && k != "altitude",
                              keys(ds))
        species = NamedTuple{Tuple(Symbol.(species_vars))}(
                      Array(ds[k]) for k in species_vars)
        return merge((; h_atm, ne, Te), species)
    end
end
