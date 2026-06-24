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
    AltitudeGrid(bottom, top; dz0=0.15, growth_scale=60, z_transition=100, dz_max=10, scale=1)

Create an altitude grid from `bottom` to `top` (in km). The keyword arguments control the
grid spacing, see [`make_altitude_grid`](@ref); in particular `scale` multiplies every step
and is meant for convergence testing.
"""
struct AltitudeGrid{FT, V<:AbstractVector{FT}} <: AbstractGrid
    h::V
    Δh::V
    n::Int
    bottom::FT
    top::FT
end

function AltitudeGrid(bottom, top; kwargs...)
    h = make_altitude_grid(bottom, top; kwargs...)
    Δh = diff(h)
    FT = eltype(h)
    return AltitudeGrid{FT, typeof(h)}(h, Δh, length(h), FT(h[1] / 1e3), FT(h[end] / 1e3))
end

"""
    make_altitude_grid(bottom_altitude, top_altitude;
                       dz0=0.15, growth_scale=60, z_transition=100, dz_max=10, scale=1)

Create an altitude grid (in m) spanning exactly `bottom_altitude` to `top_altitude` (in km).

The step size is uniform (`dz0`) below `z_transition`, then grows exponentially with
altitude up to a cap:

    dz(z) = scale * min(dz_max, dz0 * exp((z - z_transition) / growth_scale))    [km]

The fine uniform spacing at the bottom and the exponential coarsening above are matched to
where electron transport develops sharp vertical structure: the energy-deposition layers
in the lower thermosphere (~95-200 km), whose lower flanks are much steeper than the
neutral scale height. Higher up the flux profiles smooth out and coarse steps lose no
accuracy. With the defaults, an 80-700 km grid has ~580 points with steps from 150 m at
the bottom to 10 km at the top (grid-convergence measurements on steady-state runs show
this resolution carries a few % error at the peak of the volume-excitation profile, see
`scale` below).

The first and last grid points land exactly on the requested limits: the uniform step is
rounded so that `z_transition` falls on a grid point, and the small leftover at the top is
spread evenly across the graded steps (a uniform stretch of well under one step) so the
grid ends exactly at `top_altitude` without introducing an anomalously small top step.

# Arguments
- `bottom_altitude`: altitude, in km, for the bottom of the simulation
- `top_altitude`: altitude, in km, for the top of the simulation

# Keyword Arguments
- `dz0 = 0.15`: step size, in km, below `z_transition`
- `growth_scale = 60`: e-folding scale, in km, of the step size above `z_transition`
- `z_transition = 100`: altitude, in km, where the step size starts to grow
- `dz_max = 10`: maximum step size, in km (reached above ~350 km with the defaults)
- `scale = 1`: multiplies every grid step, keeping the same step-size profile in altitude.
    Use it for grid-convergence testing: rerun a case with `scale = 0.5` and compare the
    results to quantify the discretisation error of the default grid for that case. The
    transport scheme converges at first order in `scale`: halving `scale` halves the
    error, doubling `scale` doubles it.

# Outputs
- `h_atm`: altitude (m). Vector [nZ]
"""
function make_altitude_grid(bottom_altitude, top_altitude;
                            dz0 = 0.15, growth_scale = 60.0, z_transition = 100.0,
                            dz_max = 10.0, scale = 1.0)
    bottom_altitude < top_altitude || throw(ArgumentError(
        "bottom_altitude must be below top_altitude, got $(bottom_altitude) ≥ $(top_altitude)."))
    (dz0 > 0 && dz_max > 0 && scale > 0 && growth_scale > 0) || throw(ArgumentError(
        "dz0, dz_max, scale and growth_scale must all be positive."))

    # Step size (km) as a function of altitude (km)
    dz(z) = scale * min(dz_max, dz0 * exp(max(0.0, z - z_transition) / growth_scale))

    # Uniform segment, with the step rounded so that the segment ends exactly on
    # min(z_transition, top_altitude)
    knee = min(z_transition, top_altitude)
    if bottom_altitude < z_transition
        n_uniform = max(1, round(Int, (knee - bottom_altitude) / (scale * dz0)))
        h = collect(range(Float64(bottom_altitude), knee; length = n_uniform + 1))
    else
        h = [Float64(bottom_altitude)]
    end

    # Graded segment: grow steps by the law until the next one would overshoot, then land
    # exactly on top_altitude by spreading the leftover across all graded steps (a small
    # uniform stretch). This keeps the smooth step-size ratio with no special last step, so
    # the smallest step in the grid stays the bottom one (important for the CFL limit, which
    # is set by minimum(diff(z))).
    i_knee = length(h)
    while h[end] + dz(h[end]) < top_altitude
        push!(h, h[end] + dz(h[end]))
    end
    if h[end] < top_altitude
        if length(h) - i_knee ≥ 1   # at least one graded step: distribute the remainder
            z_knee = h[i_knee]
            stretch = (top_altitude - z_knee) / (h[end] - z_knee)
            for i in (i_knee + 1):length(h)
                h[i] = z_knee + (h[i] - z_knee) * stretch
            end
        else                        # top is within one step of the knee: single step up
            push!(h, Float64(top_altitude))
        end
    end

    return h .* 1e3
end

"""
    suggest_bottom_altitude(E_max, msis_file; safety=2.0)

Suggest a bottom altitude (in km) for the simulation, deep enough that electrons with
energies up to `E_max` (in eV) are fully stopped before reaching it.

The estimate combines the empirical range of electrons in air,
`R(E) = 4.30e-7 + 5.36e-6 * (E/1keV)^1.67` g cm⁻² (Rees, 1989), with the mass-density
profile of `msis_file`: the suggestion is the highest altitude at which the atmospheric
column mass depth exceeds `safety * R(E_max)`. The default `safety = 2` places the
boundary a few km below the altitude where the downward flux has decayed to ~0, so the
absorbing bottom boundary condition (Ie = 0) does not clip the solution and the full
profile of the deposition peak is captured, without wasting grid points at depth.

If the atmosphere of `msis_file` does not reach deep enough, a warning is emitted and the
lowest altitude of the file is returned; generate a deeper MSIS file to go further down.

# Arguments
- `E_max`: highest electron energy of the simulation (eV)
- `msis_file`: path to the MSIS atmosphere file

# Keyword Arguments
- `safety = 2.0`: multiplier on the electron range

# Returns
- suggested bottom altitude (km)
"""
function suggest_bottom_altitude(E_max, msis_file; safety = 2.0)
    data = load_msis_data(msis_file)
    z_km = data.height_km
    ρ = data.air                              # total mass density (kg/m³)
    # Column mass depth M(z) = ∫ᶻᵗᵒᵖ ρ dz′, converted from kg/m² to g/cm²
    M = zeros(length(z_km))
    for i in (length(z_km) - 1):-1:1
        M[i] = M[i + 1] + 0.5 * (ρ[i] + ρ[i + 1]) * (z_km[i + 1] - z_km[i]) * 1e3 * 0.1
    end
    R = electron_range_in_air(E_max)
    i_bottom = findlast(>=(safety * R), M)
    if isnothing(i_bottom)
        @warn "Electrons of $(E_max) eV can penetrate below the lowest altitude of the " *
              "MSIS file ($(z_km[1]) km). Returning $(z_km[1]) km; generate a deeper " *
              "MSIS file to place the bottom boundary lower."
        return z_km[1]
    end
    return z_km[i_bottom]
end

"""
    electron_range_in_air(E)

Empirical range of an electron with energy `E` (in eV) in air, in g cm⁻² (eq. 3.3.4 at p.40
in Chapter 3 of Rees, 1989).
"""
electron_range_in_air(E) = 4.30e-7 + 5.36e-6 * (E / 1e3)^1.67

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
