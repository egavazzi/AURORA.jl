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
    return AltitudeGrid{FT, typeof(h)}(h, Δh, length(h), FT(h[1] / 1e3), FT(h[end] / 1e3))
end

"""
    make_altitude_grid(bottom_altitude, top_altitude; dz_max=25)

Create an altitude grid based on the altitude limits given as input. It uses constant
steps of 150 m below 100 km, and a non-linear grid (steps growing with altitude, capped at
`dz_max`) above 100 km.

# Calling
`h_atm = make_altitude_grid(bottom_altitude, top_altitude; dz_max=25)`

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
    bottom_altitude < top_altitude || throw(ArgumentError(
        "bottom_altitude must be below top_altitude, got $(bottom_altitude) ≥ $(top_altitude)."))
    bottom_m = bottom_altitude * 1e3
    top_m    = top_altitude * 1e3
    z_transition = 100e3

    # Step-growth profile above the transition
    Δz(n) = min(real(150 +
                     150 / 200 * n +
                     1.2 * exp(Complex((n - 150) / 22)^0.9)),
                     dz_max * 1e3)

    # Build grid from the transition up to `top`, generating as many steps as needed
    h_atm = [z_transition]
    n = 1
    while h_atm[end] < top_m
        push!(h_atm, h_atm[end] + Δz(n))
        n += 1
    end
    h_atm = h_atm[h_atm .<= top_m]   # clip at the top

    # Bottom above transition: snap to the nearest grid point
    if bottom_m ≥ z_transition
        return h_atm[h_atm .≥ bottom_m]
    end

    # Bottom below the transition: uniform fill landing exactly on the transition
    n_sub_points = max(1, round(Int, (z_transition - bottom_m) / Δz(1)))
    sub_grid = collect(range(bottom_m, z_transition; length = n_sub_points + 1))
    return vcat(sub_grid[1:end - 1], h_atm)
end

function Base.show(io::IO, grid::AltitudeGrid)
    print(io, "AltitudeGrid($(grid.bottom) - $(grid.top) km, $(grid.n) points)")
end

function Base.show(io::IO, ::MIME"text/plain", grid::AltitudeGrid)
    println(io, "AltitudeGrid:")
    println(io, "├── Range: $(grid.bottom) - $(grid.top) km")
    println(io, "├── Points: $(grid.n)")
    println(io, "├── Min spacing: $(round(minimum(grid.Δh), digits=1)) m")
    print(io, "└── Max spacing: $(round(maximum(grid.Δh), digits=1)) m")
end
