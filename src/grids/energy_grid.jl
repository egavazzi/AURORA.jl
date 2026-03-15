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
    E_edges, ΔE = make_energy_grid(E_max)
    E_centers = E_edges[1:end-1] .+ ΔE ./ 2
    FT = eltype(E_edges)
    return EnergyGrid{FT, typeof(E_edges)}(E_edges, E_centers, ΔE, length(ΔE), FT(E_max))
end

function Base.show(io::IO, grid::EnergyGrid)
    print(io, "EnergyGrid($(grid.E_edges[1])-$(grid.E_edges[end]) eV, $(grid.n) bins)")
end

function Base.show(io::IO, ::MIME"text/plain", grid::EnergyGrid)
    println(io, "EnergyGrid:")
    println(io, "├── Range: $(round(grid.E_edges[1], digits=2)) - $(round(grid.E_edges[end], digits=2)) eV")
    println(io, "├── Bins: $(grid.n)")
    println(io, "├── Min bin width: $(round(minimum(grid.ΔE), digits=3)) eV")
    print(io, "└── Max bin width: $(round(maximum(grid.ΔE), digits=3)) eV")
end
