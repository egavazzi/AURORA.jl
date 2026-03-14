"""
    EnergyGrid{FT, V<:AbstractVector{FT}} <: AbstractGrid

Energy grid for electron transport.
"""
struct EnergyGrid{FT, V<:AbstractVector{FT}} <: AbstractGrid
    E::V
    dE::V
    n::Int
    E_max::FT
end

function EnergyGrid(E_max)
    E, dE = make_energy_grid(E_max)
    FT = eltype(E)
    return EnergyGrid{FT, typeof(E)}(E, dE, length(E), FT(E_max))
end

function Base.show(io::IO, grid::EnergyGrid)
    print(io, "EnergyGrid($(grid.E[1])-$(grid.E_max) eV, $(grid.n) bins)")
end

function Base.show(io::IO, ::MIME"text/plain", grid::EnergyGrid)
    println(io, "EnergyGrid:")
    println(io, "├── Range: $(round(grid.E[1], digits=2)) - $(grid.E_max) eV")
    println(io, "├── Bins: $(grid.n)")
    println(io, "├── Min bin width: $(round(minimum(grid.dE), digits=3)) eV")
    print(io, "└── Max bin width: $(round(maximum(grid.dE), digits=3)) eV")
end
