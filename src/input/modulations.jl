## ====================================================================================== ##
## Abstract type
## ====================================================================================== ##

"""
    AbstractModulation

Abstract supertype for all temporal modulation types.

Concrete subtypes: [`ConstantModulation`](@ref), [`SinusoidalFlickering`](@ref),
[`SquareFlickering`](@ref), [`SmoothOnset`](@ref).
"""
abstract type AbstractModulation end



## ====================================================================================== ##
## ConstantModulation
## ====================================================================================== ##

"""
    ConstantModulation()

No temporal modulation — the flux is constant in time.

Used as the default modulation for steady-state simulations and for file-based input
where the time dependence is already contained in the file.
"""
struct ConstantModulation <: AbstractModulation end

function Base.show(io::IO, ::ConstantModulation)
    print(io, "ConstantModulation()")
end



## ====================================================================================== ##
## SinusoidalFlickering
## ====================================================================================== ##

"""
    SinusoidalFlickering(f; amplitude=1.0)

Sinusoidal temporal modulation of the electron flux.

The flux oscillates between `(1 - amplitude)` and `1` of the base flux, using the
pattern `1 - cos²(πft)`. With `amplitude=1.0`, the flux modulates fully between 0 and 1.

# Arguments
- `f`: modulation frequency (Hz)

# Keyword Arguments
- `amplitude=1.0`: modulation depth. 0 = constant flux, 1 = full on/off modulation.

# Examples
```julia
SinusoidalFlickering(10.0)                  # 10 Hz, full modulation
SinusoidalFlickering(5.0; amplitude=0.5)    # 5 Hz, partial modulation
```
"""
struct SinusoidalFlickering <: AbstractModulation
    f::Float64
    amplitude::Float64

    function SinusoidalFlickering(f; amplitude=1.0)
        f > 0 || error("Frequency must be positive, got $f")
        0 <= amplitude <= 1 || error("Amplitude must be in [0, 1], got $amplitude")
        new(Float64(f), Float64(amplitude))
    end
end

function Base.show(io::IO, m::SinusoidalFlickering)
    print(io, "SinusoidalFlickering(f=$(m.f) Hz, amplitude=$(m.amplitude))")
end



## ====================================================================================== ##
## SquareFlickering
## ====================================================================================== ##

"""
    SquareFlickering(f; amplitude=1.0)

Square-wave temporal modulation of the electron flux.

The flux alternates between `(1 - amplitude)` and `1` of the base flux at the
given frequency. With `amplitude=1.0`, the flux alternates fully between 0 and 1.

# Arguments
- `f`: modulation frequency (Hz)

# Keyword Arguments
- `amplitude=1.0`: modulation depth. 0 = constant flux, 1 = full on/off modulation.

# Examples
```julia
SquareFlickering(10.0)                  # 10 Hz, full modulation
SquareFlickering(5.0; amplitude=0.5)    # 5 Hz, partial modulation
```
"""
struct SquareFlickering <: AbstractModulation
    f::Float64
    amplitude::Float64

    function SquareFlickering(f; amplitude=1.0)
        f > 0 || error("Frequency must be positive, got $f")
        0 <= amplitude <= 1 || error("Amplitude must be in [0, 1], got $amplitude")
        new(Float64(f), Float64(amplitude))
    end
end

function Base.show(io::IO, m::SquareFlickering)
    print(io, "SquareFlickering(f=$(m.f) Hz, amplitude=$(m.amplitude))")
end



## ====================================================================================== ##
## SmoothOnset
## ====================================================================================== ##

"""
    SmoothOnset(t_start, t_end)

Smooth flux onset using a C∞-smooth transition function.

The flux transitions smoothly from 0 to 1 over the interval `[t_start, t_end]`.
Before `t_start`, the flux is 0; after `t_end`, the flux is 1.

When `t_start == t_end`, this acts as a step function (Heaviside) at `t_start`.

# Arguments
- `t_start`: start time of the transition (s)
- `t_end`: end time of the transition (s)

# Examples
```julia
SmoothOnset(0.0, 0.1)    # smooth ramp from t=0 to t=0.1 s
SmoothOnset(0.0, 0.0)    # step function at t=0
```
"""
struct SmoothOnset <: AbstractModulation
    t_start::Float64
    t_end::Float64

    function SmoothOnset(t_start, t_end)
        t_end >= t_start || error("t_end ($t_end) must be >= t_start ($t_start)")
        new(Float64(t_start), Float64(t_end))
    end
end

function Base.show(io::IO, m::SmoothOnset)
    print(io, "SmoothOnset(t_start=$(m.t_start) s, t_end=$(m.t_end) s)")
end



## ====================================================================================== ##
## apply_modulation
## ====================================================================================== ##

"""
    apply_modulation(mod::AbstractModulation, t_shifted)

Apply temporal modulation to a time grid, returning a vector of modulation factors
(values between 0 and 1) for each time step.

`t_shifted` is the time grid shifted for energy- and angle-dependent electron travel times,
important in case a `z_source` was specified in the `InputFlux`.
"""
function apply_modulation end


"""
    apply_modulation(::ConstantModulation, t_shifted)

Returns ones — no temporal modulation. The flux is simply turned on at `t=0`.
"""
function apply_modulation(::ConstantModulation, t_shifted)
    return Float64.(t_shifted .>= 0)
end


"""
    apply_modulation(mod::SinusoidalFlickering, t_shifted)

Sinusoidal modulation: oscillates between `(1 - amplitude)` and `1`.
"""
function apply_modulation(mod::SinusoidalFlickering, t_shifted)
    # Base pattern: (1 - cos²(πft)) goes from 0 to 1
    base_modulation = 1.0 .- cos.(π * mod.f .* t_shifted).^2
    # Scale by amplitude and shift: result goes from (1-amplitude) to 1
    modulated = (1.0 - mod.amplitude) .+ mod.amplitude .* base_modulation
    # Apply onset (Heaviside at t=0)
    return modulated .* (t_shifted .>= 0)
end


"""
    apply_modulation(mod::SquareFlickering, t_shifted)

Square-wave modulation: alternates between `(1 - amplitude)` and `1`.
"""
function apply_modulation(mod::SquareFlickering, t_shifted)
    # Square wave helper function
    function _square(x)
        ifelse(mod2pi(x) < π, 1.0, -1.0)
    end
    # Base pattern: (1 + square(...))/2 goes from 0 to 1
    base_modulation = (1.0 .+ _square.(2π * mod.f .* t_shifted .- π/2)) ./ 2
    # Scale by amplitude and shift
    modulated = (1.0 - mod.amplitude) .+ mod.amplitude .* base_modulation
    # Apply onset (Heaviside at t=0)
    return modulated .* (t_shifted .>= 0)
end


"""
    apply_modulation(mod::SmoothOnset, t_shifted)

Smooth onset using a C∞-smooth transition function over `[t_start, t_end]`.
"""
function apply_modulation(mod::SmoothOnset, t_shifted)
    return _smooth_transition.(t_shifted, mod.t_start, mod.t_end)
end
