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

"""
    make_altitude_grid(bottom_altitude, top_altitude; dz_max=25)

Create an altitude grid based on the altitude limits given as input. It uses
constant steps of 150m for altitudes below 100km, and a non-linear grid above 100km.

# Calling
`h_atm = make_altitude_grid(bottom_altitude, top_altitude; max_dz=25)`

# Arguments
- `bottom_altitude`: altitude, in km, for the bottom of the simulation
- `top_altitude`: altitude, in km, for the top of the simulation

# Keyword Arguments
- `dz_max = 25`: maximum step size, in km. Relevant for high altitudes where the numbers
    can get large.

# Outputs
- `h_atm`: altitude (m). Vector [nZ]
"""
function make_altitude_grid(bottom_altitude, top_altitude; dz_max = 25)
    Δz(n) = 150 .+
            150 / 200 * (0:(n - 1)) .+
            1.2 * exp.(Complex.(((0:(n - 1)) .- 150) / 22) .^ 0.9)
    dz = real.(Δz(500))
    dz = min.(dz, dz_max * 1e3)
    # With the default dz_max of 25 km, the grid can go up to ~2700 km (which is way
    # too high, don't do that).
    h_atm = 100e3 .+ cumsum(dz) .- dz[1]
    h_atm = vcat((bottom_altitude * 1e3):(h_atm[2] - h_atm[1]):h_atm[1], h_atm[2:end]) # add altitude steps under 100km
    i_zmax = findlast(h_atm .<= top_altitude * 1e3)
    h_atm = h_atm[1:i_zmax]
    return h_atm
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
