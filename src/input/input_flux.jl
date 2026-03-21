using MAT: matopen
using DataInterpolations: ConstantInterpolation, LinearInterpolation


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
        if !(interpolation in (:constant, :linear))
            error("Unknown interpolation type: $interpolation. Must be :constant or :linear.")
        end

        # Resample each (beam, energy) time series onto the simulation grid.
        # Values outside the file's time range are set to zero.
        t_lo, t_hi = t_top[1], t_top[end]
        for i_μ in 1:n_μ, iE in 1:n_E
            flux_series = Float64.(Ie_top_file[i_μ, :, iE])
            if interpolation == :constant
                itp = ConstantInterpolation(flux_series, t_top)
            else  # :linear
                itp = LinearInterpolation(flux_series, t_top)
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
## Ie_top_modulated function
## ====================================================================================== ##

"""
    Ie_top_modulated(IeE_tot, model, Beams, t, n_loop;
                     spectrum=:flat, E_min=0.0, E₀=nothing, ΔE=nothing,
                     modulation=:none, f=0.0, amplitude=1.0,
                     z_source=model.altitude_grid.h[end]/1e3, t_start=0.0, t_end=0.0)

Create a time-dependent electron flux distribution with configurable energy spectrum shape
and temporal modulation.

Returns an electron flux distribution (in #e⁻/m²/s) such that when integrated over energy
(at full modulation), the total energy flux `IeE_tot` is recovered.

# Arguments
- `IeE_tot`: total energy flux (W/m²)
- `model`: `AuroraModel` describing the grids and atmosphere
- `Beams`: indices of the electron beams with precipitating flux
- `t`: time grid (s). Range or Vector [n_t]
- `n_loop`: number of loops (for repeated simulations)

# Keyword Arguments
## Energy spectrum shape
- `spectrum=:flat`: energy spectrum type. Either `:flat` or `:gaussian`
- `E_min=0.0`: minimum energy threshold (eV) - only used for `spectrum=:flat`
- `E₀=nothing`: characteristic/center energy (eV) - required for `spectrum=:gaussian`
- `ΔE=nothing`: energy width (eV) - required for `spectrum=:gaussian`

## Temporal modulation
- `modulation=:none`: temporal modulation type. One of `:none`, `:sinus`, or `:square`
- `f=0.0`: modulation frequency (Hz) - used for `:sinus` and `:square`
- `amplitude=1.0`: modulation depth. 0 = constant flux, 1 = full on/off modulation.
  Values between 0 and 1 give partial modulation.

## Time-dependent features
- `z_source`: altitude of the precipitation source (km). Default: top of ionosphere
- `t_start=0.0`: start time for smooth flux onset (s) - only used for `modulation=:none`
- `t_end=0.0`: end time for smooth flux onset (s) - only used for `modulation=:none`

# Returns
- `Ie_top`: electron number flux (#e⁻/m²/s). Matrix [n_beams, n_t, n_E]

# Physics
The function creates an electron precipitation spectrum with:

1. **Energy spectrum**: Either flat (uniform in #e⁻/m²/s/eV above E_min) or Gaussian
   (centered at E₀ with width ΔE). Both are normalized so the total energy flux equals
   `IeE_tot`.

2. **Energy- and angle-dependent arrival times**: Lower energy electrons travel slower,
   arriving later at the ionosphere top. Similarly, electrons with high pitch-angles
   travel longer paths. This creates natural dispersion.

3. **Temporal modulation**: The flux can be constant (`:none`), sinusoidally modulated
   (`:sinus`), or square-wave modulated (`:square`). The `amplitude` parameter controls
   the modulation depth.

# Examples
Flat spectrum with smooth onset:
```jldoctest
julia> model = AuroraModel((100, 600), 180:-10:0, 10e3, find_msis_file(), find_iri_file());

julia> Ie = Ie_top_modulated(1e-2, model, 1:2, 0:0.01:1, 1;
                             spectrum=:flat, E_min=9000.0, t_start=0.0, t_end=0.1);
```

Gaussian spectrum with sinusoidal modulation at 10 Hz:
```jldoctest
julia> msis_file = find_msis_file(verbose=false);

julia> iri_file = find_iri_file(verbose=false);

julia> model = AuroraModel((100, 600), 180:-10:0, 10e3, msis_file, iri_file; verbose=false);

julia> Ie = Ie_top_modulated(1e-2, model, 1:2, 0:0.001:0.5, 1;
                             spectrum=:gaussian, E₀=5000.0, ΔE=500.0,
                             modulation=:sinus, f=10.0, amplitude=1.0);
```
"""
function Ie_top_modulated(IeE_tot, model::AuroraModel, Beams, t, n_loop;
                          spectrum=:flat, E_min=0.0, E₀=nothing, ΔE=nothing,
                          modulation=:none, f=0.0, amplitude=1.0,
                          z_source=model.altitude_grid.h[end]/1e3, t_start=0.0, t_end=0.0)

    E_centers = model.energy_grid.E_centers
    ΔE_grid = model.energy_grid.ΔE
    μ_center = model.pitch_angle_grid.μ_center
    Ω_beam = model.scattering.Ω_beam
    z = model.altitude_grid.h

    ## ==================== Input validation ==================== ##
    if spectrum == :flat
        if E_min < 0
            error("E_min is negative ($E_min eV). E_min must be a positive value.")
        end
        if E_min > E_centers[end]
            error("E_min ($E_min eV) is larger than the maximum energy in the grid ($(E_centers[end]) eV).")
        end
    elseif spectrum == :gaussian
        if isnothing(E₀) || isnothing(ΔE)
            error("For spectrum=:gaussian, both E₀ and ΔE must be specified.")
        end
        if E₀ <= 0 || ΔE <= 0
            error("E₀ ($E₀ eV) and ΔE ($ΔE eV) must be positive values.")
        end
    else
        error("Unknown spectrum type: $spectrum. Must be :flat or :gaussian.")
    end

    if !(modulation in (:none, :sinus, :square))
        error("Unknown modulation type: $modulation. Must be :none, :sinus, or :square.")
    end

    if z_source < z[end] / 1e3
        error("z_source ($z_source km) is below the top of the simulated ionosphere \
        ($(z[end]/1e3) km). It must be above or equal.")
    end

    if amplitude < 0 || amplitude > 1
        error("The amplitude ($amplitude) is outside the range [0, 1].")
    end

    ## ==================== Setup ==================== ##
    # Physical constants
    qₑ = 1.602176620898e-19  # Elementary charge (C)

    # Calculate total number of time steps
    n_time = (n_loop - 1) * (length(t) - 1) + length(t)

    # Initialize output array [n_beams, n_time, n_energy]
    Ie_top = zeros(length(μ_center), n_time, length(E_centers))

    # Convert total energy flux from W/m² to eV/m²/s
    IeE_tot_eV = IeE_tot / qₑ

    # Compute energy spectrum (differential flux per eV)
    if spectrum == :flat
        Φ_spectrum = _flat_spectrum(IeE_tot_eV, E_centers, ΔE_grid, E_min)
    else  # :gaussian
        Φ_spectrum = _gaussian_spectrum(IeE_tot_eV, E_centers, ΔE_grid, E₀, ΔE)
    end

    # Create extended time grid for all loops
    t_total = t[1]:Float64(t.step):(t[end] * n_loop)

    # Calculate distance from source to top of ionosphere (m)
    z_distance = z_source * 1e3 - z[end]

    # Calculate reference time shift (for highest energy in first beam)
    # Used to align all arrival times relative to the fastest electrons
    t_ref = z_distance / (abs(μ_center[Beams[1]]) * v_of_E(E_centers[end]))

    ## ==================== Main loop ==================== ##
    for i_μ in Beams
        if μ_center[i_μ] >= 0
            continue  # Skip upward-going beams
        end
        beam_fraction = Ω_beam[i_μ] / sum(Ω_beam[Beams])

        for iE in eachindex(E_centers)
            # Skip if no flux at this energy
            if Φ_spectrum[iE] ≈ 0
                continue
            end

            # Calculate base flux for this beam and energy bin
            flux_base = Φ_spectrum[iE] * beam_fraction * ΔE_grid[iE]

            # Calculate travel time for electrons of this energy and pitch angle
            t_travel = z_distance / (abs(μ_center[i_μ]) * v_of_E(E_centers[iE]))

            # Time-shifted grid: subtract travel time difference relative to reference
            t_shifted = t_total .- (t_travel - t_ref)

            # Apply temporal modulation
            temporal_factor = _apply_modulation(t_shifted, modulation, f, amplitude,
                                                t_start, t_end)

            # Assign flux
            Ie_top[i_μ, :, iE] = flux_base .* temporal_factor
        end
    end

    return Ie_top
end


"""
    Ie_top_modulated(IeE_tot, model, Beams;
                     spectrum=:flat, E_min=0.0, E₀=nothing, ΔE=nothing)

Steady-state version of `Ie_top_modulated` that does not include time-dependent behavior.

This overload is designed for steady-state simulations and eliminates the need to specify
time-related parameters (`t`, `n_loop`, `modulation`, `f`, `amplitude`, `z_source`, etc.).
It internally calls the time-dependent version with minimal time grid (`1:1:1`) and `n_loop=1`.

# Arguments
- `IeE_tot`: total energy flux (W/m²)
- `model`: `AuroraModel` describing the grids and atmosphere
- `Beams`: indices of the electron beams with precipitating flux

# Keyword Arguments
- `spectrum=:flat`: energy spectrum type. Either `:flat` or `:gaussian`
- `E_min=0.0`: minimum energy threshold (eV) - only used for `spectrum=:flat`
- `E₀=nothing`: characteristic/center energy (eV) - required for `spectrum=:gaussian`
- `ΔE=nothing`: energy width (eV) - required for `spectrum=:gaussian`

# Returns
- `Ie_top`: electron number flux (#e⁻/m²/s). Matrix [n_beams, 1, n_E]

# Examples
```jldoctest
julia> msis_file = find_msis_file(verbose=false);

julia> iri_file = find_iri_file(verbose=false);

julia> model = AuroraModel((100, 600), 180:-10:0, 10e3, msis_file, iri_file; verbose=false);

julia> Ie = Ie_top_modulated(1e-2, model, 1:2; spectrum=:flat, E_min=9000.0);
```
"""
function Ie_top_modulated(IeE_tot, model::AuroraModel, Beams;
                          spectrum=:flat, E_min=0.0, E₀=nothing, ΔE=nothing)
    z = model.altitude_grid.h
    # Call time-dependent version with dummy time grid and default time-related parameters
    return Ie_top_modulated(IeE_tot, model, Beams,
                            1:1:1, 1;
                           spectrum=spectrum, E_min=E_min, E₀=E₀, ΔE=ΔE,
                           modulation=:none, f=0.0, amplitude=1.0,
                           z_source=z[end]/1e3, t_start=0.0, t_end=0.0)
end


## ====================================================================================== ##
## Ie_with_LET function – for steady-state only
## ====================================================================================== ##

"""
    Ie_with_LET(IeE_tot, E₀, model, Beams; low_energy_tail=true)

Return an electron spectra following a Maxwellian distribution with a low
energy tail (LET)

This function is a **corrected** implementation of Meier/Strickland/Hecht/Christensen
JGR 1989 (pages 13541-13552)

# Arguments
- `IeE_tot`: total energy flux (W/m²)
- `E₀`: characteristic energy (eV)
- `model`: `AuroraModel` describing the grids and atmosphere
- `Beams`: indices of the electron beams with a precipitating flux

# Keyword Arguments
- `low_energy_tail=true`: control the presence of a low energy tail

# Returns:
- `Ie_top`: electron number flux (#e⁻/m²/s). Matrix [n_beams, 1, nE]

# Important notes
This is a corrected version of the equations present in Meier et al. 1989
to match the results presented in Fig. 4 of their paper.\\
Changes were made to the factor `b`:
- no inverse

# Examples:
Calling the function with flux only in the two first beams (0 to 20°) and an "isotropic"
pitch-angle distribution.
```jldoctest
julia> model = AuroraModel((100, 600), 180:-10:0, 10e3, find_msis_file(), find_iri_file());

julia> Ie = AURORA.Ie_with_LET(1e-2, 1e3, model, 1:2);

```

Calling the function with flux only in the three first beams (0 to 30°) and a
custom pitch-angle distribution (1/2 of the total flux in the first beam,
1/4 in the second beam and 1/4 in the third beam).
```jldoctest
julia> msis_file = find_msis_file(verbose=false);

julia> iri_file = find_iri_file(verbose=false);

julia> model = AuroraModel((100, 600), 180:-10:0, 10e3, msis_file, iri_file; verbose=false);

julia> Ie = Ie_with_LET(1e-2, 1e3, model, 1:3);
```
"""
function Ie_with_LET(IeE_tot, E₀, model::AuroraModel, Beams; low_energy_tail=true)
    E_centers = model.energy_grid.E_centers
    ΔE = model.energy_grid.ΔE
    μ_center = model.pitch_angle_grid.μ_center
    Ω_beam = model.scattering.Ω_beam

    Ie_top = zeros(length(μ_center), 1, length(E_centers))

    # Physical constants
    qₑ = 1.602176620898e-19  # Elementary charge (C)

    # Convert total energy flux from W/m² to eV/m²/s
    IeE_tot_eV = IeE_tot / qₑ  # eV/m²/s

    # Maxwellian spectra
    # π is gone as we do not normalize in /ster
    # However we keep the 2 as it seems necessary to keep the total energy flux
    # Either the Meier et al. conversion to ster is slightly wrong, or there is something I don't get?
    Φₘ = IeE_tot_eV / (2 * E₀^3) .* E_centers .* exp.(-E_centers ./ E₀)
    # Parameter for the LET (corrected equations to match the Fig. 4)
    b = (0.8 * E₀) .* (E₀ < 500) +
        (0.1 * E₀ + 350) .* (E₀ >= 500) # use 350 instead of .35 because we are in eV

    for i_μ in Beams
        if low_energy_tail
            Ie_max = maximum(Φₘ) # max of the Maxwellian - to scale LET amplitude
            Ie_top[i_μ, 1, :] = (Φₘ .+ 0.4 * Ie_max * (E₀ ./ E_centers) .* exp.(-E_centers ./ b)) .*
                                ΔE * Ω_beam[i_μ] ./ sum(Ω_beam[Beams])
        else
            Ie_top[i_μ, 1, :] = Φₘ .*
                                ΔE * Ω_beam[i_μ] ./ sum(Ω_beam[Beams])
        end
    end

    return Ie_top
end



## ====================================================================================== ##
## Helper functions for Ie_top_modulated
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


"""
    _apply_modulation(t_shifted, modulation, f, amplitude, t_start, t_end)

Apply temporal modulation to the flux based on the shifted time grid.

# Arguments
- `t_shifted`: time grid shifted for energy/angle-dependent delays
- `modulation`: `:none`, `:sinus`, or `:square`
- `f`: frequency for sinus/square modulation (Hz)
- `amplitude`: modulation depth (0 = constant, 1 = full on/off)
- `t_start`, `t_end`: smooth onset interval (only used for `:none`)

# Returns
- Vector of modulation factors (0 to 1) for each time step
"""
function _apply_modulation(t_shifted, modulation, f, amplitude, t_start, t_end)
    if modulation == :none
        # Smooth onset/offset using smooth_transition
        return _smooth_transition.(t_shifted, t_start, t_end)
    elseif modulation == :sinus
        # Sinusoidal modulation: oscillates between (1-amplitude) and 1
        # Base pattern: (1 - cos²(πft)) goes from 0 to 1
        base_modulation = 1.0 .- cos.(π * f .* t_shifted).^2
        # Scale by amplitude and shift: result goes from (1-amplitude) to 1
        modulated = (1.0 - amplitude) .+ amplitude .* base_modulation
        # Apply onset (Heaviside at t=0)
        return modulated .* (t_shifted .>= 0)
    elseif modulation == :square
        # Square wave modulation: oscillates between (1-amplitude) and 1
        # Base pattern: (1 + square(...))/2 goes from 0 to 1
        base_modulation = (1.0 .+ square.(2π * f .* t_shifted .- π/2)) ./ 2
        # Scale by amplitude and shift
        modulated = (1.0 - amplitude) .+ amplitude .* base_modulation
        # Apply onset (Heaviside at t=0)
        return modulated .* (t_shifted .>= 0)
    else
        error("Unknown modulation type: $modulation. Must be :none, :sinus, or :square.")
    end
end

function square(x)
    ifelse(mod2pi(x) < π, 1.0, -1.0)
end
