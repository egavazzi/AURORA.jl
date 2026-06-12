using NCDatasets: NCDataset

"""
    VolumeExcitationResult

Result of the volume excitation rate computation. Contains the volume excitation/ionization
rates on the altitude x time grid, as saved in `analysis/volume_excitation.nc`.

# Fields
- `Q4278`, `Q6730`, `Q7774`, `Q8446`, `QO1D`, `QO1S`: volume excitation rates (photons/m³/s) as `Matrix{Float64}` (n_z x n_t)
- `Q7774_O`, `Q7774_O2`: species-resolved contributions to `Q7774`
- `Q8446_O`, `Q8446_O2`: species-resolved contributions to `Q8446`
- `QOi`, `QO2i`, `QN2i`: ionization rates (/m³/s) as `Matrix{Float64}` (n_z x n_t)
- `h_atm`: altitude grid (m) as `Vector{Float64}`
- `t`: time grid (s) as `Vector{Float64}`
- `savedir`: directory from which the result was loaded, or `nothing`
"""
struct VolumeExcitationResult
    Q4278::Matrix{Float64}
    Q6730::Matrix{Float64}
    Q7774::Matrix{Float64}
    Q7774_O::Matrix{Float64}
    Q7774_O2::Matrix{Float64}
    Q8446::Matrix{Float64}
    Q8446_O::Matrix{Float64}
    Q8446_O2::Matrix{Float64}
    QO1D::Matrix{Float64}
    QO1S::Matrix{Float64}
    QOi::Matrix{Float64}
    QO2i::Matrix{Float64}
    QN2i::Matrix{Float64}
    h_atm::Vector{Float64}
    t::Vector{Float64}
    savedir::Union{String, Nothing}
end

"""
    ColumnExcitationResult

Result of the column-integrated excitation rate computation. Contains the column-integrated
emission intensities as a function of time, as saved in `analysis/column_excitation.nc`.

# Fields
- `I_4278`, `I_6730`, `I_7774`, `I_8446`, `I_O1D`, `I_O1S`: column-integrated intensities (photons/m²/s) as `Vector{Float64}` (n_t)
- `I_7774_O`, `I_7774_O2`: species-resolved contributions to `I_7774`
- `I_8446_O`, `I_8446_O2`: species-resolved contributions to `I_8446`
- `t`: time grid (s) as `Vector{Float64}`
"""
struct ColumnExcitationResult
    I_4278::Vector{Float64}
    I_6730::Vector{Float64}
    I_7774::Vector{Float64}
    I_7774_O::Vector{Float64}
    I_7774_O2::Vector{Float64}
    I_8446::Vector{Float64}
    I_8446_O::Vector{Float64}
    I_8446_O2::Vector{Float64}
    I_O1D::Vector{Float64}
    I_O1S::Vector{Float64}
    t::Vector{Float64}
end

"""
    load_volume_excitation(directory::String)

Load the volume excitation rates from `analysis/volume_excitation.nc` in the given directory.
Returns a [`VolumeExcitationResult`](@ref).
"""
function load_volume_excitation(directory::String)
    filepath = joinpath(directory, "analysis", "volume_excitation.nc")
    NCDataset(filepath, "r") do ds
        return VolumeExcitationResult(
            Matrix{Float64}(ds["Q4278"][:, :]),
            Matrix{Float64}(ds["Q6730"][:, :]),
            Matrix{Float64}(ds["Q7774"][:, :]),
            Matrix{Float64}(ds["Q7774_O"][:, :]),
            Matrix{Float64}(ds["Q7774_O2"][:, :]),
            Matrix{Float64}(ds["Q8446"][:, :]),
            Matrix{Float64}(ds["Q8446_O"][:, :]),
            Matrix{Float64}(ds["Q8446_O2"][:, :]),
            Matrix{Float64}(ds["QO1D"][:, :]),
            Matrix{Float64}(ds["QO1S"][:, :]),
            Matrix{Float64}(ds["QOi"][:, :]),
            Matrix{Float64}(ds["QO2i"][:, :]),
            Matrix{Float64}(ds["QN2i"][:, :]),
            Vector{Float64}(ds["altitude"][:]),
            Vector{Float64}(ds["time"][:]),
            directory,
        )
    end
end

"""
    load_volume_excitation(sim::AuroraSimulation)

Load the volume excitation rates from the simulation's save directory.
"""
load_volume_excitation(sim::AuroraSimulation) = load_volume_excitation(sim.output.savedir)

"""
    load_column_excitation(directory::String)

Load the column-integrated excitation intensities from `analysis/column_excitation.nc`
in the given directory. Returns a [`ColumnExcitationResult`](@ref).
"""
function load_column_excitation(directory::String)
    filepath = joinpath(directory, "analysis", "column_excitation.nc")
    NCDataset(filepath, "r") do ds
        return ColumnExcitationResult(
            Vector{Float64}(ds["I_4278"][:]),
            Vector{Float64}(ds["I_6730"][:]),
            Vector{Float64}(ds["I_7774"][:]),
            Vector{Float64}(ds["I_7774_O"][:]),
            Vector{Float64}(ds["I_7774_O2"][:]),
            Vector{Float64}(ds["I_8446"][:]),
            Vector{Float64}(ds["I_8446_O"][:]),
            Vector{Float64}(ds["I_8446_O2"][:]),
            Vector{Float64}(ds["I_O1D"][:]),
            Vector{Float64}(ds["I_O1S"][:]),
            Vector{Float64}(ds["time"][:]),
        )
    end
end

"""
    load_column_excitation(sim::AuroraSimulation)

Load the column-integrated excitation intensities from the simulation's save directory.
"""
load_column_excitation(sim::AuroraSimulation) = load_column_excitation(sim.output.savedir)

"""
    IeTopResult

Electron flux at the top altitude of the model, as saved in `analysis/Ie_top.nc` by
[`make_Ie_top_file`](@ref). It contains all beams — both downward (precipitation) and
upward (backscattered) — and is derived from the simulation *output*, so it is distinct
from the `Ie_input` boundary condition stored in `simulation_data.nc`.

# Fields
- `Ietop`: flux array of size (n_beams x n_t x n_E), units: m⁻² s⁻¹
- `t`: time grid (s), relative so always start at 0
- `E_edges`: energy bin edges (eV), length n_E + 1
- `E_centers`: energy bin centres (eV), length n_E
- `ΔE`: energy bin widths (eV), length n_E
- `μ_lims`: cosine-of-pitch-angle limits, length n_beams + 1
"""
struct IeTopResult
    Ietop::Array{Float64, 3}
    t::Vector{Float64}
    E_edges::Vector{Float64}
    E_centers::Vector{Float64}
    ΔE::Vector{Float64}
    μ_lims::Vector{Float64}
end

"""
    load_Ie_top(directory::String)

Load the electron flux at the top altitude of the model from `analysis/Ie_top.nc`
in `directory`. Contains all beams (downward precipitation and upward backscatter).
Returns an [`IeTopResult`](@ref).
"""
function load_Ie_top(directory::String)
    filepath = joinpath(directory, "analysis", "Ie_top.nc")
    NCDataset(filepath, "r") do ds
        t = Vector{Float64}(ds["time"][:])
        t = t .- t[1]   # make relative (start at 0)
        E_edges   = Vector{Float64}(ds["energy_edges"][:])
        E_centers = Vector{Float64}(ds["energy"][:])
        ΔE = diff(E_edges)
        Ietop = Array{Float64, 3}(ds["Ie_top_raw"][:, :, :])
        μ_lims = Vector{Float64}(ds["mu_lims"][:])
        return IeTopResult(Ietop, t, E_edges, E_centers, ΔE, μ_lims)
    end
end
