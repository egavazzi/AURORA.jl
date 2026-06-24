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
steps below 100 km, and a non-linear grid (steps growing with altitude, capped at
`dz_max`) above 100 km.

# Calling
`h_atm = make_altitude_grid(bottom_altitude, top_altitude; dz_max=25)`

# Arguments
- `bottom_altitude`: altitude, in km, for the bottom of the simulation
- `top_altitude`: altitude, in km, for the top of the simulation

# Keyword Arguments
- `dz_max = 25`: maximum step size, in km. Relevant for high altitudes where the numbers
    can get large.

# Notes
- A `bottom_altitude` above 100 km is honoured (the grid snaps to the nearest grid point
    at or above it) instead of silently collapsing back to 100 km.
- The uniform segment below 100 km lands exactly on the transition, so there is no
    anomalous step straddling 100 km.
- The grid extends as far up as needed to reach `top_altitude`; there is no hard-coded
    cap on the number of points.

# Outputs
- `h_atm`: altitude (m). Vector [nZ]
"""
function make_altitude_grid(bottom_altitude, top_altitude; dz_max = 25)
    bottom_altitude < top_altitude || throw(ArgumentError(
        "bottom_altitude must be below top_altitude, got $(bottom_altitude) ≥ $(top_altitude)."))
    bot_m, top_m, z_transition = bottom_altitude * 1e3, top_altitude * 1e3, 100e3

    # Legacy step-growth profile above the transition; step number `n` = 1, 2, …, capped at
    # `dz_max`. The shape is unchanged from the previous version — only the assembly around
    # it is fixed below.
    step_at(n) = min(real(150 + 150 / 200 * n + 1.2 * exp(Complex((n - 150) / 22)^0.9)),
                     dz_max * 1e3)

    # Graded segment from the transition up to `top`, generating as many steps as needed
    # (the old code hard-coded 500 steps, which capped the point count and the reachable
    # altitude).
    graded = [z_transition]
    n = 1
    while graded[end] < top_m
        push!(graded, graded[end] + step_at(n))
        n += 1
    end
    graded = graded[graded .<= top_m]   # clip at the top (legacy behaviour, no resizing)

    if bot_m ≥ z_transition
        # Bottom in the F-region: snap to the nearest grid point at or above it, instead of
        # silently collapsing the grid back to 100 km.
        return graded[graded .≥ bot_m]
    end

    # Bottom below the transition: uniform fill landing exactly on the transition (or on
    # `top`, if that is itself below 100 km), removing the old spacing anomaly at 100 km.
    knee = min(z_transition, top_m)
    n_sub = max(1, round(Int, (knee - bot_m) / step_at(1)))
    sub = collect(range(bot_m, knee; length = n_sub + 1))
    return knee < top_m ? vcat(sub[1:end - 1], graded) : sub
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
