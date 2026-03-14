"""
    AltitudeGrid{FT, V<:AbstractVector{FT}} <: AbstractGrid

Altitude grid for the ionosphere simulation.

# Fields
- `h`: altitude values (m)
- `Δh`: grid spacing (m)
- `n`: number of grid points
- `bottom`: bottom altitude (km)
- `top`: top altitude (km)

# Constructor
    AltitudeGrid(bottom, top; dz_max=25)

Create an altitude grid from `bottom` to `top` (in km), with maximum grid spacing `dz_max` (km).
"""
struct AltitudeGrid{FT, V<:AbstractVector{FT}} <: AbstractGrid
    h::V
    Δh::V
    n::Int
    bottom::FT
    top::FT
end

function AltitudeGrid(bottom, top; dz_max=25)
    h = make_altitude_grid(bottom, top; dz_max=dz_max)
    Δh = diff(h)
    FT = eltype(h)
    return AltitudeGrid{FT, typeof(h)}(h, Δh, length(h), FT(bottom), FT(top))
end

function Base.show(io::IO, grid::AltitudeGrid)
    print(io, "AltitudeGrid($(grid.bottom)-$(grid.top) km, $(grid.n) points)")
end

function Base.show(io::IO, ::MIME"text/plain", grid::AltitudeGrid)
    println(io, "AltitudeGrid:")
    println(io, "├── Range: $(grid.bottom) - $(grid.top) km")
    println(io, "├── Points: $(grid.n)")
    println(io, "├── Min spacing: $(round(minimum(grid.Δh), digits=1)) m")
    print(io, "└── Max spacing: $(round(maximum(grid.Δh), digits=1)) m")
end
