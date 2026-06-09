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
    load_coordinates(sim_dir) → NamedTuple

Load everything in [`SimulationResult`](@ref) *except* the bulky `Ie` array: the coordinate
vectors (`t`, `h_atm`, `E_centers`, `E_edges`, `ΔE`, `μ_lims`, `beam_weights`), the `Ie`
dimensions `(n_z, n_μ, n_t, n_E)`, and `savedir`. Cheap to call.
"""
function load_coordinates(sim_dir::AbstractString)
    nc_path = joinpath(sim_dir, "simulation_data.nc")
    NCDataset(nc_path, "r") do ds
        n_z, n_μ, n_t, n_E = size(ds["Ie"])
        E_edges = Vector{Float64}(ds["energy_edges"][:])
        return (; t            = Vector{Float64}(ds["time"][:]),
                  h_atm        = Vector{Float64}(ds["altitude"][:]),
                  E_centers    = Vector{Float64}(ds["energy"][:]),
                  E_edges,
                  ΔE           = diff(E_edges),
                  μ_lims       = Vector{Float64}(ds["mu_lims"][:]),
                  beam_weights = Vector{Float64}(ds["beam_weight"][:]),
                  n_z, n_μ, n_t, n_E,
                  savedir      = sim_dir)
    end
end


"""
    foreach_Ie_time_chunk(f, sim_dir; max_bytes = 512 * 1024^2)

Stream the electron flux `Ie` from `simulation_data.nc` in `sim_dir` over contiguous
time-chunks, calling `f(Ie_chunk, t_range)` for each. `Ie_chunk::Array{Float64,4}` has shape
`[n_z, n_μ, length(t_range), n_E]` and `t_range` is the corresponding range of time indices.

The chunk length is the largest that keeps one `Ie_chunk` under `max_bytes`, so peak memory is
bounded regardless of the total run length. On disk `Ie` is chunked one time-slice per chunk,
so each slab read is contiguous.
"""
function foreach_Ie_time_chunk(f, sim_dir::AbstractString; max_bytes::Real = 512 * 1024^2)
    nc_path = joinpath(sim_dir, "simulation_data.nc")
    NCDataset(nc_path, "r") do ds
        n_z, n_μ, n_t, n_E = size(ds["Ie"])
        slice_bytes = n_z * n_μ * n_E * sizeof(Float64)
        nt_chunk = time_chunk_length(slice_bytes, max_bytes, n_t)
        for t0 in 1:nt_chunk:n_t
            t_range = t0:min(t0 + nt_chunk - 1, n_t)
            f(Array{Float64, 4}(ds["Ie"][:, :, t_range, :]), t_range)
        end
    end
    return nothing
end


"""
    load_results(sim_dir::AbstractString; tidx=:, zidx=:, μidx=:, eidx=:, max_bytes=2*1024^3)
    load_results(sim::AuroraSimulation; kwargs...)

Load the raw electron-flux output from `simulation_data.nc` in `sim_dir` (or in the
simulation's save directory). Returns a [`SimulationResult`](@ref).

The whole `Ie` array is read into memory eagerly. For large time-dependent runs this can be
many GB; use the `tidx`/`zidx`/`μidx`/`eidx` selectors to load only a sub-box, or stream the
data with [`foreach_Ie_time_chunk`](@ref). As a safety net `load_results` refuses to allocate
more than `max_bytes` (default 2 GiB) and reports how to narrow the load; pass `max_bytes=Inf`
to disable the guard.

# Selectors
Each selector is a `Colon` (`:`, the default — load everything) or a contiguous integer range
(e.g. `tidx = 1:100`). The matching coordinate vectors and edge-derived quantities (`E_edges`,
`ΔE`, `μ_lims`) are sliced consistently. Non-contiguous indexing is unsupported; read the
`NCDataset` directly if you need it.

# Example
```julia
res = load_results("my_run")              # everything
res = load_results("my_run"; tidx=1:100)  # first 100 time slices only
res.Ie                                     # flux [n_z, n_μ, n_t, n_E]
res.t                                      # time axis (s)
```
"""
function load_results(sim_dir::AbstractString;
                      tidx = Colon(), zidx = Colon(), μidx = Colon(), eidx = Colon(),
                      max_bytes::Real = 2 * 1024^3)
    nc_path = joinpath(sim_dir, "simulation_data.nc")
    NCDataset(nc_path, "r") do ds
        n_z, n_μ, n_t, n_E = size(ds["Ie"])
        zsel = resolve_selector(zidx, n_z, "zidx")
        μsel = resolve_selector(μidx, n_μ, "μidx")
        tsel = resolve_selector(tidx, n_t, "tidx")
        esel = resolve_selector(eidx, n_E, "eidx")

        nbytes = length(zsel) * length(μsel) * length(tsel) * length(esel) * sizeof(Float64)
        if nbytes > max_bytes
            throw(ArgumentError(
                "load_results would allocate $(round(nbytes / 1024^3; digits=2)) GiB " *
                "(> max_bytes = $(round(max_bytes / 1024^3; digits=2)) GiB). Narrow the load " *
                "with the tidx/zidx/μidx/eidx selectors, raise max_bytes (or set it to Inf), " *
                "or stream the data with foreach_Ie_time_chunk."))
        end

        Ie           = Array{Float64, 4}(ds["Ie"][zsel, μsel, tsel, esel])  # [n_z, n_μ, n_t, n_E]
        t            = Vector{Float64}(ds["time"][tsel])
        h_atm        = Vector{Float64}(ds["altitude"][zsel])
        E_centers    = Vector{Float64}(ds["energy"][esel])
        E_edges      = bounds_for(Vector{Float64}(ds["energy_edges"][:]), esel)
        ΔE           = diff(E_edges)
        μ_lims       = bounds_for(Vector{Float64}(ds["mu_lims"][:]), μsel)
        beam_weights = Vector{Float64}(ds["beam_weight"][μsel])
        return SimulationResult(Ie, t, h_atm, E_centers, E_edges, ΔE, μ_lims, beam_weights, sim_dir)
    end
end

load_results(sim::AuroraSimulation; kwargs...) = load_results(sim.output.savedir; kwargs...)


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


# Resolve a user selector (`:` or a contiguous integer range) into a concrete range,
# validating bounds. Non-contiguous indexing is rejected so the edge-derived quantities
# (`E_edges`, `μ_lims`, `ΔE`) stay well-defined; drop to the `NCDataset` directly for that.
function resolve_selector(sel, n::Integer, name::AbstractString)
    sel === Colon() && return 1:n
    if sel isa AbstractUnitRange{<:Integer}
        (!isempty(sel) && first(sel) >= 1 && last(sel) <= n) ||
            throw(ArgumentError("$name = $sel is out of bounds for dimension of length $n"))
        return sel
    end
    throw(ArgumentError(
        "$name must be a Colon (:) or a contiguous integer range (e.g. 1:10); got " *
        "$(typeof(sel)). For non-contiguous indexing, read the NCDataset directly."))
end

# Slice an (n+1)-length bounds vector (energy edges, μ limits) to match a contiguous
# selection of the n-length centered dimension.
bounds_for(bounds, sel::AbstractUnitRange) = bounds[first(sel):(last(sel) + 1)]

# Number of time slices per streaming chunk: the largest count whose `slice_bytes` stays
# under `max_bytes`, clamped to `[1, n_t]`. Handles `max_bytes = Inf` (→ all `n_t` at once).
function time_chunk_length(slice_bytes::Real, max_bytes::Real, n_t::Integer)
    per_chunk = max_bytes / slice_bytes
    return per_chunk >= n_t ? Int(n_t) : max(1, floor(Int, per_chunk))
end
