"""
    EnergyGrid{FT, V<:AbstractVector{FT}} <: AbstractGrid

Energy grid for electron transport.

# Fields
- `E_edges`: energy bin edges (eV), length `n + 1`
- `E_centers`: energy bin centers (eV), length `n`
- `ΔE`: energy bin widths (eV), length `n`
- `n`: number of energy bins
- `E_max`: requested maximum energy (eV)
"""
struct EnergyGrid{FT, V<:AbstractVector{FT}} <: AbstractGrid
    E_edges::V
    E_centers::V
    ΔE::V
    n::Int
    E_max::FT
end

function EnergyGrid(E_max)
    E_edges, E_centers, ΔE = make_energy_grid(E_max)
    FT = eltype(E_edges)
    return EnergyGrid{FT, typeof(E_edges)}(E_edges, E_centers, ΔE, length(ΔE), FT(E_max))
end

"""
    make_energy_grid(E_max)

Create an energy grid based on the maximum energy `E_max` given as input.

# Calling
`E_edges, E_centers, ΔE = make_energy_grid(E_max)`

# Inputs
- `E_max`: upper limit for the energy grid (in eV)

# Outputs
- `E_edges`: energy bin edges (eV). Vector [nE + 1]
- `E_centers`: energy bin centers (eV). Vector [nE]
- `ΔE`: energy bin widths (eV). Vector [nE]
"""
function make_energy_grid(E_max)
    if E_max > 1e6
        error("AURORA does not support energies above 1 MeV. Please use a lower value for E_max.")
    end
    E_function(X, dE_initial, dE_final, C, X0) = dE_initial + (1 + tanh(C * (X - X0))) / 2 * dE_final
    E = cumsum(E_function.(0:100000, 0.15, 11.5, 0.05, 80)) .+ 1.9
    iE_max = findmin(abs.(E .- E_max))[2]  # find the index for the upper limit of the energy grid
    E = E[1:iE_max]                        # crop E accordingly
    ΔE = diff(E)
    ΔE = [ΔE; ΔE[end]]
    E_edges = [E; E[end] + ΔE[end]]         # add the last upper edge
    E_centers = E_edges[1:end-1] .+ ΔE ./ 2
    return E_edges, E_centers, ΔE
end

function Base.show(io::IO, grid::EnergyGrid)
    print(io, "EnergyGrid($(round(grid.E_edges[1], digits=2)) - $(round(grid.E_edges[end], digits=2)) eV, $(grid.n) bins)")
end

function Base.show(io::IO, ::MIME"text/plain", grid::EnergyGrid)
    println(io, "EnergyGrid:")
    println(io, "├── Range: $(round(grid.E_edges[1], digits=2)) - $(round(grid.E_edges[end], digits=2)) eV")
    println(io, "├── Bins: $(grid.n)")
    println(io, "├── Min bin width: $(round(minimum(grid.ΔE), digits=3)) eV")
    print(io, "└── Max bin width: $(round(maximum(grid.ΔE), digits=3)) eV")
end
