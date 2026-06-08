using NCDatasets: NCDataset

# ======================================================================================== #
#                              RESULT FILE READERS                                       #
# ======================================================================================== #

"""
    SimulationResult

The raw electron-flux output of a simulation, as loaded from `simulation_data.nc` by
[`load_results`](@ref).

# Fields
- `Ie`        : electron number flux `[n_z, n_μ, n_t, n_E]` (m⁻² s⁻¹). See the note on
                [Output & data](@ref Output): this is a *number* flux, already integrated over
                each energy bin and beam solid angle.
- `t`         : time axis `[n_t]` (s)
- `h_atm`     : altitude grid `[n_z]` (m)
- `E_centers` : energy bin centres `[n_E]` (eV)
- `E_edges`   : energy bin edges `[n_E+1]` (eV)
- `ΔE`          : energy bin widths `[n_E]` (eV) (= `diff(E_edges)`)
- `μ_lims`      : pitch-angle cosine bin boundaries `[n_μ+1]`
- `beam_weights`: solid-angle beam weights `[n_μ]` (sr)
- `savedir`     : directory the result was loaded from, or `nothing`
"""
struct SimulationResult
    Ie::Array{Float64, 4}
    t::Vector{Float64}
    h_atm::Vector{Float64}
    E_centers::Vector{Float64}
    E_edges::Vector{Float64}
    ΔE::Vector{Float64}
    μ_lims::Vector{Float64}
    beam_weights::Vector{Float64}
    savedir::Union{String, Nothing}
end

function Base.show(io::IO, r::SimulationResult)
    n_z, n_μ, n_t, n_E = size(r.Ie)
    print(io, "SimulationResult(Ie ", n_z, "×", n_μ, "×", n_t, "×", n_E,
          ", ", n_t, " time steps)")
end

"""
    load_results(sim_dir::AbstractString) → SimulationResult
    load_results(sim::AuroraSimulation)   → SimulationResult

Load the raw electron-flux output from `simulation_data.nc` in `sim_dir` (or in the
simulation's save directory). Returns a [`SimulationResult`](@ref).

# Example
```julia
res = load_results("my_run")
res.Ie     # flux [n_z, n_μ, n_t, n_E]
res.t      # time axis (s)
```
"""
function load_results(sim_dir::AbstractString)
    nc_path = joinpath(sim_dir, "simulation_data.nc")
    NCDataset(nc_path, "r") do ds
        Ie        = Array{Float64, 4}(ds["Ie"][:, :, :, :])   # [n_z, n_μ, n_t, n_E]
        t         = Vector{Float64}(ds["time"][:])
        h_atm     = Vector{Float64}(ds["altitude"][:])
        E_centers = Vector{Float64}(ds["energy"][:])
        E_edges      = Vector{Float64}(ds["energy_edges"][:])
        ΔE           = diff(E_edges)
        μ_lims       = Vector{Float64}(ds["mu_lims"][:])
        beam_weights = Vector{Float64}(ds["beam_weight"][:])
        return SimulationResult(Ie, t, h_atm, E_centers, E_edges, ΔE, μ_lims, beam_weights, sim_dir)
    end
end

load_results(sim::AuroraSimulation) = load_results(sim.output.savedir)


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
