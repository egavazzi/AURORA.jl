using MAT: matread

"""
    VolumeExcitationResult

Result of the volume excitation rate computation. Contains the volume excitation/ionization
rates on the altitude x time grid, as saved in `Qzt_all_L.mat`.

# Fields
- `Q4278`, `Q6730`, `Q7774`, `Q8446`, `QO1D`, `QO1S`: volume excitation rates (photons/m³/s) as `Matrix{Float64}` (n_z × n_t)
- `Q7774_O`, `Q7774_O2`: species-resolved contributions to `Q7774`
- `Q8446_O`, `Q8446_O2`: species-resolved contributions to `Q8446`
- `QOi`, `QO2i`, `QN2i`: ionization rates (/m³/s) as `Matrix{Float64}` (n_z × n_t)
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
emission intensities as a function of time, as saved in `I_lambda_of_t.mat`.

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

Load the volume excitation rates from `Qzt_all_L.mat` in the given directory.
Returns a [`VolumeExcitationResult`](@ref).
"""
function load_volume_excitation(directory::String)
    filepath = joinpath(directory, "Qzt_all_L.mat")
    data = matread(filepath)
    return VolumeExcitationResult(
        data["Q4278"],
        data["Q6730"],
        data["Q7774"],
        data["Q7774_O"],
        data["Q7774_O2"],
        data["Q8446"],
        data["Q8446_O"],
        data["Q8446_O2"],
        data["QO1D"],
        data["QO1S"],
        data["QOi"],
        data["QO2i"],
        data["QN2i"],
        vec(collect(data["h_atm"])),
        vec(collect(data["t"])),
        directory,
    )
end

"""
    load_volume_excitation(sim::AuroraSimulation)

Load the volume excitation rates from the simulation's save directory.
"""
load_volume_excitation(sim::AuroraSimulation) = load_volume_excitation(sim.savedir)

"""
    load_column_excitation(directory::String)

Load the column-integrated excitation intensities from `I_lambda_of_t.mat` in the given directory.
Returns a [`ColumnExcitationResult`](@ref).
"""
function load_column_excitation(directory::String)
    filepath = joinpath(directory, "I_lambda_of_t.mat")
    data = matread(filepath)
    return ColumnExcitationResult(
        vec(collect(data["I_4278"])),
        vec(collect(data["I_6730"])),
        vec(collect(data["I_7774"])),
        vec(collect(data["I_7774_O"])),
        vec(collect(data["I_7774_O2"])),
        vec(collect(data["I_8446"])),
        vec(collect(data["I_8446_O"])),
        vec(collect(data["I_8446_O2"])),
        vec(collect(data["I_O1D"])),
        vec(collect(data["I_O1S"])),
        vec(collect(data["t"])),
    )
end

"""
    load_column_excitation(sim::AuroraSimulation)

Load the column-integrated excitation intensities from the simulation's save directory.
"""
load_column_excitation(sim::AuroraSimulation) = load_column_excitation(sim.savedir)

"""
    IeTopResult

Incoming electron flux at the top of the ionosphere, as saved in `Ie_incoming_*.mat`.

# Fields
- `Ietop`: flux array of size (n_beams x n_t x n_E)
- `t`: time grid (s), relative so always start at 0
- `E_edges`: energy bin edges (eV), length n_E + 1
- `E_centers`: energy bin centres (eV), length n_E
- `ΔE`: energy bin widths (eV), length n_E
- `mu_lims`: cosine-of-pitch-angle limits, length n_beams + 1
"""
struct IeTopResult
    Ietop::Array{Float64, 3}
    t::Vector{Float64}
    E_edges::Vector{Float64}
    E_centers::Vector{Float64}
    ΔE::Vector{Float64}
    mu_lims::Vector{Float64}
end

"""
    load_input(directory::String)

Load the incoming electron flux at the top of the ionosphere from `Ie_incoming_*.mat`
in `directory`. Returns an [`IeTopResult`](@ref).
"""
function load_input(directory::String)
    filepath = find_input_file(directory)
    data = matread(filepath)
    t = vec(collect(data["t_top"]))
    t = [t; t[end] + diff(t)[end]] .- t[1]
    E_edges =  vec(collect(data["E_edges"]))
    E_centers = vec(collect(data["E_centers"]))
    ΔE = vec(collect(data["dE"]))
    return IeTopResult(
        data["Ie_total"],
        t,
        E_edges,
        E_centers,
        ΔE,
        vec(collect(data["mu_lims"])),
    )
end
