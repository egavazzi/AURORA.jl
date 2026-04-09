using MAT: matopen
using DataInterpolations: ConstantInterpolation, LinearInterpolation


## ====================================================================================== ##
## InputFlux struct
## ====================================================================================== ##

"""
    InputFlux{S<:AbstractSpectrum, M<:AbstractModulation}

A unified electron precipitation input, combining an energy spectrum shape
with a temporal modulation and beam selection.

# Fields
- `spectrum`: the energy spectrum shape (any [`AbstractSpectrum`](@ref))
- `modulation`: the temporal modulation (any [`AbstractModulation`](@ref))
- `beams`: indices of the electron beams with precipitating flux
- `z_source`: altitude of the precipitation source (km)

# Constructors
```julia
# Steady-state (constant modulation is the default)
InputFlux(FlatSpectrum(1e-2; E_min=2900); beams=1:2)
InputFlux(MaxwellianSpectrum(1e-2, 1000); beams=1:2)

# Time-dependent with modulation
InputFlux(FlatSpectrum(1e-2; E_min=100), SinusoidalFlickering(5.0);
          beams=1, z_source=3000)

# From file (modulation must be and is by default ConstantModulation)
InputFlux(FileSpectrum("path.mat"))
```
"""
struct InputFlux{S<:AbstractSpectrum, M<:AbstractModulation}
    spectrum::S
    modulation::M
    beams::Vector{Int}
    z_source::Float64
end

# Convenience constructor: spectrum only (constant modulation)
function InputFlux(spectrum::AbstractSpectrum; beams=1, z_source=NaN)
    InputFlux(spectrum, ConstantModulation(); beams=beams, z_source=z_source)
end

# Main keyword constructor: spectrum + modulation
function InputFlux(spectrum::AbstractSpectrum, modulation::AbstractModulation;
                   beams=1, z_source=NaN)
    if spectrum isa FileSpectrum && !(modulation isa ConstantModulation)
        error("FileSpectrum does not support temporal modulation. " *
              "The file already contains the complete flux. " *
              "Use ConstantModulation() (the default).")
    end
    beams_vec = _to_beam_vector(beams)
    InputFlux{typeof(spectrum), typeof(modulation)}(
        spectrum, modulation, beams_vec, Float64(z_source))
end

# Internal: convert any beam specification to Vector{Int}
_to_beam_vector(b::Int) = [b]
_to_beam_vector(b::AbstractRange{<:Integer}) = collect(Int, b)
_to_beam_vector(b::AbstractVector{<:Integer}) = convert(Vector{Int}, b)


function Base.show(io::IO, flux::InputFlux)
    print(io, "InputFlux($(flux.spectrum), $(flux.modulation), beams=$(flux.beams))")
end

function Base.show(io::IO, ::MIME"text/plain", flux::InputFlux)
    println(io, "InputFlux:")
    println(io, "├── Spectrum:   ", flux.spectrum)
    println(io, "├── Modulation: ", flux.modulation)
    println(io, "├── Beams:      ", flux.beams)
    z_str = isnan(flux.z_source) ? "top of ionosphere" : "$(flux.z_source) km"
    print(io,   "└── Source altitude: ", z_str)
end



## ====================================================================================== ##
## Ie_top_from_file (used by compute_flux with FileSpectrum)
## ====================================================================================== ##

"""
    Ie_top_from_file(t, E, μ_center, filename; interpolation=:constant)

Load a time-dependent electron flux from a `.mat` file and resample it onto the
simulation time grid `t`.

# Arguments
- `t`: full simulation time grid (s). Range or Vector [n_t]
- `E`: energy grid (eV). Vector [n_E]
- `μ_center`: electron beams average pitch angle cosine. Vector [n_μ]
- `filename`: path to the `.mat` file containing the flux data

# Keyword Arguments
- `interpolation=`: interpolation scheme used to resample the file's time grid
  onto the simulation time grid. Either:
  - `:constant` (default): each file value is held constant until the next sample.
  - `:linear`: linear interpolation between consecutive file samples.
  - `:pchip`: piecewise cubic Hermite interpolation.

# Returns
- `Ie_top`: electron number flux (#e⁻/m²/s). Array [n_μ, n_t, n_E]

# File format
The `.mat` file must contain:
- `Ie_total`: flux array of shape `[n_μ, n_t_file, n_E]`
- `t_top`: time grid (s) of the file, as a vector `[n_t_file]`. Optional when `n_t_file == 1`.

# Time grid handling
The file's time grid (`t_top`) does not need to match the simulation time grid.
It is treated as **relative**: `t_top[1]` is always aligned with `t[1]`, regardless of the
absolute values stored in the file. Only the **duration** and **spacing** of `t_top` matter.
Any combination of time step sizes and grid lengths is supported:
- **Different dt**: the file is resampled onto the simulation grid via interpolation.
- **File finer than simulation**: the flux is point-sampled, which may miss rapid variations.
- **File shorter than simulation**: flux is set to zero for times beyond the file duration.
- **File longer than simulation**: the file is simply evaluated up to `t[end]`.

# Notes
- Upward-going beams (`μ_center > 0`) are always set to zero.
- The energy dimension of the file may be larger than `length(E)`; only the first
  `length(E)` bins are used. It cannot be smaller.
"""
function Ie_top_from_file(t, E_centers, μ_center, filename; interpolation=:constant)
    ## Load file
    file = matopen(filename)
        Ie_top_file = read(file, "Ie_total")
        # if t_top not present, assume steady-state (constant flux)
        t_top = haskey(file, "t_top") ? vec(read(file, "t_top")) : [0.0]
    close(file)

    ## Check dimensions
    if size(Ie_top_file, 1) != length(μ_center)
        error("""
        The incoming flux `Ie_total` has the wrong number of beams (dim 1). \
        Got $(size(Ie_top_file, 1)), expected $(length(μ_center)). \
        Remember that Ie_total should have the shape [n_μ, n_t, n_E].""")
    end
    if size(Ie_top_file, 3) < length(E_centers)
        error("""
        The incoming flux `Ie_total` has too few energy bins (dim 3). \
        Got $(size(Ie_top_file, 3)), expected at least $(length(E_centers)). \
        Remember that Ie_total should have the shape [n_μ, n_t, n_E].""")
    end
    if length(t_top) > 1 && size(Ie_top_file, 2) != length(t_top)
        error("""
        The time dimension of `Ie_total` (dim 2, $(size(Ie_top_file, 2)) steps) does not \
        match the length of `t_top` ($(length(t_top)) steps).""")
    end

    ## Resample onto simulation grid
    n_μ = length(μ_center)
    n_t = length(t)
    n_E = length(E_centers)
    Ie_top = zeros(n_μ, n_t, n_E)

    if length(t_top) == 1
        Ie_top .= Ie_top_file[:, 1:1, 1:n_E] # this repeats the unique time slice from Ie_top_file
    else
        # Shift t_top so it starts at t[1], i.e treat file times as relative to simulation start
        t_top = t_top .- t_top[1] .+ t[1]

        t_duration_file = t_top[end] - t_top[1]
        t_duration_sim  = t[end]    - t[1]
        if t_duration_file < t_duration_sim
            @warn """The time grid in the file covers $(t_duration_file) s, which is shorter \
            than the simulation duration $(t_duration_sim) s. \
            The precipitating flux will be set to zero for t > t[1] + $(t_duration_file) s."""
        end

        dt_file = t_top[2] - t_top[1]
        dt_sim  = t[2] - t[1]
        if dt_file < dt_sim
            @warn """The file time resolution (dt=$(dt_file) s) is finer than the \
            simulation time resolution (dt=$(dt_sim) s). The flux will be point-sampled \
            onto the simulation grid, which may miss rapid variations in the file."""
        end

        # Validate interpolation keyword
        if !(interpolation in (:constant, :linear, :pchip))
            error("Unknown interpolation type: $interpolation. Must be :constant, :linear, or :pchip.")
        end

        # Resample each (beam, energy) time series onto the simulation grid.
        # Values outside the file's time range are set to zero.
        t_lo, t_hi = t_top[1], t_top[end]
        for i_μ in 1:n_μ, iE in 1:n_E
            flux_series = Float64.(Ie_top_file[i_μ, :, iE])
            if interpolation == :constant
                itp = ConstantInterpolation(flux_series, t_top)
            elseif interpolation == :linear
                itp = LinearInterpolation(flux_series, t_top)
            elseif interpolation == :pchip
                itp = PCHIPInterpolation(flux_series, t_top)
            end
            for (k, tk) in enumerate(t)
                Ie_top[i_μ, k, iE] = (tk < t_lo || tk > t_hi) ? 0.0 : itp(tk)
            end
        end
    end

    # Set the input up-flux to zero
    Ie_top[μ_center .> 0, :, :] .= 0

    return Ie_top
end



## ====================================================================================== ##
## Helper functions for spectrum evaluation
## ====================================================================================== ##

"""
    _flat_spectrum(IeE_tot_eV, E, dE, E_min)

Compute a flat (constant) differential number flux spectrum above E_min, normalized so that
the total energy flux equals IeE_tot_eV.

Returns a vector of differential number flux per eV (#e⁻/m²/s/eV) for each energy bin.
The flux is zero below E_min.
"""
function _flat_spectrum(IeE_tot_eV, E_centers, ΔE, E_min)
    # Find first index where bin lower edge <= E_min
    i_Emin = findlast((E_centers .- ΔE./2) .<= E_min)
    if isnothing(i_Emin)
        i_Emin = 1
    end

    # Calculate normalization: ∑ E_centers[i] * ΔE[i] for i >= i_Emin
    energy_integral = sum(E_centers[i_Emin:end] .* ΔE[i_Emin:end])

    # Differential number flux (constant above E_min)
    Φ = zeros(length(E_centers))
    Φ[i_Emin:end] .= IeE_tot_eV / energy_integral

    return Φ  # #e⁻/m²/s/eV
end


"""
    _gaussian_spectrum(IeE_tot_eV, E, dE, E₀, ΔE)

Compute a Gaussian differential number flux spectrum centered at E₀ with width ΔE,
normalized so that the total energy flux equals IeE_tot_eV.

The Gaussian shape is: Φ(E) ∝ exp(-(E - E₀)² / ΔE²)

Returns a vector of differential number flux per eV (#e⁻/m²/s/eV) for each energy bin.
"""
function _gaussian_spectrum(IeE_tot_eV, E_centers, ΔE, E₀, ΔE_gauss)
    # Unnormalized Gaussian shape
    Φ_shape = exp.(-(E_centers .- E₀).^2 ./ ΔE_gauss^2)

    # Calculate normalization: we want ∑ Φ[i] * E_centers[i] * ΔE[i] = IeE_tot_eV
    # So: normalization = IeE_tot_eV / ∑(Φ_shape[i] * E_centers[i] * ΔE[i])
    energy_integral = sum(Φ_shape .* E_centers .* ΔE)

    # Normalized differential number flux
    Φ = Φ_shape .* (IeE_tot_eV / energy_integral)

    return Φ  # #e⁻/m²/s/eV
end

## ====================================================================================== ##
## compute_flux — unified interface
## ====================================================================================== ##

"""
    compute_flux(flux::InputFlux, model::AuroraModel, t)

Compute the full 3D electron flux array for a **time-dependent** simulation.

Evaluates the energy spectrum, distributes it over beams, and applies temporal
modulation (including energy- and angle-dependent electron travel-time delays if a
`z_source` was specified in the InputFlux).

# Arguments
- `flux`: an [`InputFlux`](@ref) describing the precipitation
- `model`: an [`AuroraModel`](@ref) with grids and atmosphere/ionosphere
- `t`: simulation time grid (s). Range or Vector

# Returns
- `Ie_top`: electron number flux (#e⁻/m²/s). Array `[n_beams, n_t, n_E]`
"""
function compute_flux(flux::InputFlux{<:FileSpectrum}, model::AuroraModel, t)
    spec = flux.spectrum
    return Ie_top_from_file(t, model.energy_grid.E_centers,
                            model.pitch_angle_grid.μ_center, spec.filename;
                            interpolation=spec.interpolation)
end

function compute_flux(flux::InputFlux{<:AbstractSpectrum}, model::AuroraModel, t)
    E_centers = model.energy_grid.E_centers
    ΔE = model.energy_grid.ΔE
    μ_center = model.pitch_angle_grid.μ_center
    Ω_beam = model.scattering.Ω_beam
    z = model.altitude_grid.h

    # Resolve z_source: use top of ionosphere if NaN
    z_source = isnan(flux.z_source) ? z[end] / 1e3 : flux.z_source

    if z_source < z[end] / 1e3
        error("z_source ($z_source km) is below the top of the simulated ionosphere " *
              "($(z[end]/1e3) km). It must be above or equal.")
    end

    ## ==================== Setup ==================== ##
    # Initialize output array [n_beams, n_t, n_energy]
    Ie_top = zeros(length(μ_center), length(t), length(E_centers))

    # Evaluate the energy spectrum shape
    Φ_spectrum = evaluate_spectrum(flux.spectrum, model)

    # Calculate distance from source to top of ionosphere (m)
    z_distance = z_source * 1e3 - z[end]

    # Calculate reference time shift (for highest energy in first beam)
    t_ref = z_distance / (abs(μ_center[flux.beams[1]]) * v_of_E(E_centers[end]))

    ## ==================== Main loop ==================== ##
    for i_μ in flux.beams
        if μ_center[i_μ] >= 0
            continue  # Skip upward-going beams
        end
        beam_fraction = Ω_beam[i_μ] / sum(Ω_beam[flux.beams])

        for iE in eachindex(E_centers)
            if Φ_spectrum[iE] ≈ 0
                continue
            end

            # Base flux for this beam and energy bin
            flux_base = Φ_spectrum[iE] * beam_fraction * ΔE[iE]

            # Travel time for electrons of this energy and pitch angle
            t_travel = z_distance / (abs(μ_center[i_μ]) * v_of_E(E_centers[iE]))

            # Time-shifted grid: subtract travel time difference relative to reference
            t_shifted = t .- (t_travel - t_ref)

            # Apply temporal modulation
            temporal_factor = apply_modulation(flux.modulation, t_shifted)

            # Assign flux
            Ie_top[i_μ, :, iE] = flux_base .* temporal_factor
        end
    end

    return Ie_top
end


"""
    compute_flux(flux::InputFlux, model::AuroraModel)

Compute the electron flux array for a **steady-state** simulation.

This method is for steady-state runs: it evaluates the energy spectrum and distributes
it over beams without any time dependence or travel-time delays.

Only [`ConstantModulation`](@ref) is allowed for steady-state. Other modulation types
will raise an error.

# Arguments
- `flux`: an [`InputFlux`](@ref) describing the precipitation
- `model`: an [`AuroraModel`](@ref) with grids and atmosphere

# Returns
- `Ie_top`: electron number flux (#e⁻/m²/s). Array `[n_beams, 1, n_E]`
"""
function compute_flux(flux::InputFlux{<:FileSpectrum}, model::AuroraModel)
    spec = flux.spectrum
    return Ie_top_from_file(1:1:1, model.energy_grid.E_centers,
                            model.pitch_angle_grid.μ_center, spec.filename;
                            interpolation=spec.interpolation)
end

function compute_flux(flux::InputFlux{<:AbstractSpectrum, ConstantModulation}, model::AuroraModel)
    E_centers = model.energy_grid.E_centers
    ΔE = model.energy_grid.ΔE
    μ_center = model.pitch_angle_grid.μ_center
    Ω_beam = model.scattering.Ω_beam

    # Evaluate the energy spectrum shape
    Φ_spectrum = evaluate_spectrum(flux.spectrum, model)

    # Initialize output array [n_beams, 1, n_energy]
    Ie_top = zeros(length(μ_center), 1, length(E_centers))

    for i_μ in flux.beams
        if μ_center[i_μ] >= 0
            continue  # Skip upward-going beams
        end
        beam_fraction = Ω_beam[i_μ] / sum(Ω_beam[flux.beams])
        Ie_top[i_μ, 1, :] = Φ_spectrum .* beam_fraction .* ΔE
    end

    return Ie_top
end

function compute_flux(flux::InputFlux{<:AbstractSpectrum, <:AbstractModulation}, model::AuroraModel)
    error("Steady-state compute_flux does not support $(typeof(flux.modulation)). " *
          "Use ConstantModulation() for steady-state simulations, " *
          "or provide a time grid for time-dependent simulations: compute_flux(flux, model, t).")
end
