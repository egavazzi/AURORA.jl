using MAT: matopen

function Ie_top_from_old_matlab_file(t, E, n_loop, μ_center, filename)
    Ie_top = Array{Float64}(undef, length(μ_center), (n_loop - 1) * (length(t) - 1) + length(t), length(E))

    file = matopen(filename)
    Ie_top_raw = read(file, "Ie_total")
    close(file)

    # for constant input flux (e.g. first run), we need to resize the matrix from
    # [1, n_μ * [n_E, 1]] to [1, n_μ * [n_E, n_t]]
    if size(Ie_top_raw[1], 2) == 1
        for i_μ in eachindex(μ_center)
            Ie_top_raw[i_μ] = Ie_top_raw[i_μ] * ones(1, length(t) + (n_loop - 1) * (length(t) - 1)) # (e-/m²/s)
        end
    end

    # then we resize the matrix from [1, n_μ * [n_E, n_t]] to [n_μ, n_t, n_E] to be consistent with
    # the other flux matrices
    for i_μ in eachindex(μ_center)
        Ie_top[i_μ, :,  :] = Ie_top_raw[i_μ][1:length(E), :]' # (e-/m²/s)
    end

    return Ie_top
end


function Ie_top_from_ketchup(t, E, n_loop, μ_center, filename)
    Ie_top = Array{Float64}(undef, length(μ_center), (n_loop - 1) * (length(t) - 1) + length(t), length(E))

    file = matopen(filename)
    Ie_top_raw = read(file, "Ie_total")
    close(file)

    # for constant input flux (e.g. first run), we need to resize the matrix from
    # [1, n_μ * [n_E, 1]] to [1, n_μ * [n_E, n_t]]
    if size(Ie_top_raw[1], 2) == 1 && length(t) > 1
        for i_μ in eachindex(μ_center)
            Ie_top_raw[i_μ] = repeat(Ie_top_raw[i_μ], outer=(1, length(t) + (n_loop - 1) * (length(t) - 1)))
        end
    end

    # then we resize the matrix from [1, n_μ * [n_E, n_t]] to [n_μ, n_t, n_E] to be consistent with
    # the other flux matrices
    for i_μ in 1:Int(length(μ_center) / 2) # down-flux
        Ie_top[i_μ, :,  :] = Ie_top_raw[i_μ][1:length(E), :]' # (e-/m²/s)
    end
    # set the input up-flux to 0
    Ie_top[μ_center .> 0, :, :] .= 0



    return Ie_top
end



function Ie_top_from_file(t, E, μ_center, n_loop, filename)
    Nt = (n_loop - 1) * (length(t) - 1) + length(t)

    ## load the file
    file = matopen(filename)
        Ie_top_raw = read(file, "Ie_total")
        t_top = read(file, "t_top")
    close(file)

    ## check that Ie_top is matching our simulation grid
    if size(Ie_top_raw, 1) != length(μ_center)
        error("""The incoming flux Ie_top is wrongly dimensioned. Check the θ dimension. \
        Remember that Ie_top should have the shape [n_μ, n_t, n_E].""")
    end
    if size(Ie_top_raw, 3) < length(E)
        error("""The incoming flux Ie_top is wrongly dimensioned. Check the E dimension. \
        Remember that Ie_top should have the shape [n_μ, n_t, n_E].""")
        error_grid = 1
    end

    ## check the time grid of Ie_top
    if length(t_top) == 1
        # we assume constant input flux and resize the matrix from
        # [n_μ, 1, n_E] to [n_μ, n_t, n_E]
        Ie_top = repeat(Ie_top_raw, outer=(1, Nt, 1))[:, :, 1:length(E)]
    else
        dt_top = t_top[2] - t_top[1]
        dt = t[2] - t[1]
        if dt_top > dt
            # the resolution in Ie_top is coarser than our simulation, we need to repeat elements
            dt_factor = dt_top / dt
            if !(dt_factor ≈ round(dt_factor))
                error("""Problem with the time resolution. The ratio of the dt of the \
                incoming flux over the dt of the simulation is not an integer. \n
                dt_incoming = $dt_top \n
                dt_simulation = $dt \n
                ratio = $dt_factor""")
            end
            dt_factor = round(Int, dt_factor)
            Ie_top_raw = repeat(Ie_top_raw, inner=(1, dt_factor, 1))
        elseif dt_top < dt
            # the resolution in Ie_top is finer than our simulation, we need to drop elements
            dt_factor = dt / dt_top
            if !(dt_factor ≈ round(dt_factor))
                error("""Problem with the time resolution. The ratio of the dt of the \
                simulation over the dt of the incoming flux is not an integer. \n
                dt_incoming = $dt_top \n
                dt_simulation = $dt \n
                ratio = $dt_factor""")
            end
            dt_factor = round(Int, dt_factor)
            Ie_top_raw = Ie_top_raw[:, 1:dt_factor:end, :]
        end

        if size(Ie_top_raw, 2) > Nt
            # in that case the array is too long and we need to cut it
            Ie_top_raw = Ie_top_raw[:, 1:Nt, :]
        elseif size(Ie_top_raw, 2) < Nt
            # in that case the array is too short and we fill it with zeros
            missing_data = zeros(length(μ_center), Nt - size(Ie_top_raw, 2), size(Ie_top_raw, 3))
            Ie_top_raw = cat(Ie_top_raw, missing_data; dims=2)
        end
        Ie_top = Ie_top_raw[:, :, 1:length(E)]
    end

    # set the input up-flux to 0
    Ie_top[μ_center .> 0, :, :] .= 0

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
function _flat_spectrum(IeE_tot_eV, E, dE, E_min)
    E_mid = E .+ dE ./ 2

    # Find first index where E >= E_min
    i_Emin = findlast(E .<= E_min)
    if isnothing(i_Emin)
        i_Emin = 1
    end

    # Calculate normalization: ∑ E_mid[i] * dE[i] for i >= i_Emin
    energy_integral = sum(E_mid[i_Emin:end] .* dE[i_Emin:end])

    # Differential number flux (constant above E_min)
    Φ = zeros(length(E))
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
function _gaussian_spectrum(IeE_tot_eV, E, dE, E₀, ΔE)
    E_mid = E .+ dE ./ 2

    # Unnormalized Gaussian shape
    Φ_shape = exp.(-(E_mid .- E₀).^2 ./ ΔE^2)

    # Calculate normalization: we want ∑ Φ[i] * E_mid[i] * dE[i] = IeE_tot_eV
    # So: normalization = IeE_tot_eV / ∑(Φ_shape[i] * E_mid[i] * dE[i])
    energy_integral = sum(Φ_shape .* E_mid .* dE)

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


## ====================================================================================== ##
## Main Ie_top_modulated function
## ====================================================================================== ##

"""
    Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, BeamWeight, t, n_loop, h_atm;
                     spectrum=:flat, E_min=0.0, E₀=nothing, ΔE=nothing,
                     modulation=:none, f=0.0, amplitude=1.0,
                     z_source=h_atm[end]/1e3, t_start=0.0, t_end=0.0)

Create a time-dependent electron flux distribution with configurable energy spectrum shape
and temporal modulation.

Returns an electron flux distribution (in #e⁻/m²/s) such that when integrated over energy
(at full modulation), the total energy flux `IeE_tot` is recovered.

# Arguments
- `IeE_tot`: total energy flux (W/m²)
- `E`: energy grid (eV). Vector [nE]
- `dE`: energy bin widths (eV). Vector [nE]
- `μ_center`: electron beams average pitch angle cosine. Vector [n_beams]
- `Beams`: indices of the electron beams with precipitating flux
- `BeamWeight`: weights of the different beams. Vector [n_beams]
- `t`: time grid (s). Range or Vector [n_t]
- `n_loop`: number of loops (for repeated simulations)
- `h_atm`: altitude grid (m). Vector [n_z]

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
- `Ie_top`: differential electron number flux (#e⁻/m²/s). Matrix [n_beams, n_t, n_E]

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
julia> E, dE = make_energy_grid(10e3);

julia> θ_lims = 180:-10:0;

julia> μ_center = mu_avg(θ_lims);

julia> BeamWeight = beam_weight(θ_lims);

julia> h_atm = make_altitude_grid(100, 600);

julia> Ie = Ie_top_modulated(1e-2, E, dE, μ_center, 1:2, BeamWeight, 0:0.01:1, 1, h_atm;
                             spectrum=:flat, E_min=9000.0, t_start=0.0, t_end=0.1);
```

Gaussian spectrum with sinusoidal modulation at 10 Hz:
```jldoctest
julia> E, dE = make_energy_grid(10e3);

julia> θ_lims = 180:-10:0;

julia> μ_center = mu_avg(θ_lims);

julia> BeamWeight = beam_weight(θ_lims);

julia> h_atm = make_altitude_grid(100, 600);

julia> Ie = Ie_top_modulated(1e-2, E, dE, μ_center, 1:2, BeamWeight, 0:0.001:0.5, 1, h_atm;
                             spectrum=:gaussian, E₀=5000.0, ΔE=500.0,
                             modulation=:sinus, f=10.0, amplitude=1.0);
```
"""
function Ie_top_modulated(IeE_tot, E, dE, μ_center, Beams, BeamWeight, t, n_loop, h_atm;
                          spectrum=:flat, E_min=0.0, E₀=nothing, ΔE=nothing,
                          modulation=:none, f=0.0, amplitude=1.0,
                          z_source=h_atm[end]/1e3, t_start=0.0, t_end=0.0)

    ## ==================== Input validation ==================== ##
    if spectrum == :flat
        if E_min < 0
            error("E_min is negative ($E_min eV). E_min must be a positive value.")
        end
        if E_min > E[end] + dE[end]/2
            error("E_min ($E_min eV) is larger than the maximum energy in the grid ($(E[end]) eV).")
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

    if z_source < h_atm[end] / 1e3
        error("z_source ($z_source km) is below the top of the simulated ionosphere \
        ($(h_atm[end]/1e3) km). It must be above or equal.")
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
    Ie_top = zeros(length(μ_center), n_time, length(E))

    # Convert total energy flux from W/m² to eV/m²/s
    IeE_tot_eV = IeE_tot / qₑ

    # Compute energy spectrum (differential flux per eV)
    if spectrum == :flat
        Φ_spectrum = _flat_spectrum(IeE_tot_eV, E, dE, E_min)
    else  # :gaussian
        Φ_spectrum = _gaussian_spectrum(IeE_tot_eV, E, dE, E₀, ΔE)
    end

    # Create extended time grid for all loops
    t_total = t[1]:Float64(t.step):(t[end] * n_loop)

    # Calculate distance from source to top of ionosphere (m)
    z_distance = z_source * 1e3 - h_atm[end]

    # Calculate reference time shift (for highest energy in first beam)
    # Used to align all arrival times relative to the fastest electrons
    t_ref = z_distance / (abs(μ_center[Beams[1]]) * v_of_E(E[end]))

    ## ==================== Main loop ==================== ##
    for i_μ in Beams
        if μ_center[i_μ] >= 0
            continue  # Skip upward-going beams
        end
        beam_fraction = BeamWeight[i_μ] / sum(BeamWeight[Beams])

        for iE in eachindex(E)
            # Skip if no flux at this energy
            if Φ_spectrum[iE] ≈ 0
                continue
            end

            # Calculate base flux for this beam and energy bin
            flux_base = Φ_spectrum[iE] * beam_fraction * dE[iE]

            # Calculate travel time for electrons of this energy and pitch angle
            t_travel = z_distance / (abs(μ_center[i_μ]) * v_of_E(E[iE]))

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
    Ie_with_LET(IeE_tot, E₀, E, dE, μ_center, BeamWeight, Beams; low_energy_tail=true)

Return an electron spectra following a Maxwellian distribution with a low
energy tail (LET)

This function is a **corrected** implementation of Meier/Strickland/Hecht/Christensen
JGR 1989 (pages 13541-13552)

# Arguments
- `IeE_tot`: total energy flux (W/m²)
- `E₀`: characteristic energy (eV)
- `E`: energy grid (eV). Vector [nE]
- `dE`: energy bin sizes(eV). Vector [nE]
- `μ_center`: electron beams average pitch angle cosine. Vector [n_beams]
- `BeamWeight`: weights of the different beams. Vector [n_beams]
- `Beams`: indices of the electron beams with a precipitating flux
- `low_energy_tail=true`: control the presence of a low energy tail

# Returns:
- `Ie_top`: differential electron number flux (#e⁻/m²/s). Matrix [n_beams, 1, nE]

# Important notes
This is a corrected version of the equations present in Meier et al. 1989
to match the results presented in Fig. 4 of their paper.\\
Changes were made to the factor `b`:
- no inverse

# Examples:
Calling the function with flux only in the two first beams (0 to 20°) and an "isotropic"
pitch-angle distribution.
```jldoctest
julia> E, dE = make_energy_grid(100e3);

julia> θ_lims = 180:-10:0;

julia> μ_center = mu_avg(θ_lims);

julia> BeamWeight = beam_weight(180:-10:0);

julia> Ie = AURORA.Ie_with_LET(1e-2, 1e3, E, dE, μ_center, BeamWeight, 1:2);

```

Calling the function with flux only in the three first beams (0 to 30°) and a
custom pitch-angle distribution (1/2 of the total flux in the first beam,
1/4 in the second beam and 1/4 in the third beam).
```jldoctest
julia> E, dE = make_energy_grid(100e3);

julia> θ_lims = 180:-10:0;

julia> μ_center = mu_avg(θ_lims);

julia> BeamWeight = [2, 1, 1];

julia> Ie = Ie_with_LET(1e-2, 1e3, E, dE, μ_center, BeamWeight, 1:3);

```
"""
function Ie_with_LET(IeE_tot, E₀, E, dE, μ_center, BeamWeight, Beams; low_energy_tail=true)
    Ie_top = zeros(length(μ_center), 1, length(E))

    # Physical constants
    qₑ = 1.602176620898e-19  # Elementary charge (C)

    # Convert total energy flux from W/m² to eV/m²/s
    IeE_tot_eV = IeE_tot / qₑ  # eV/m²/s

    E_middle = E .+ dE / 2 # middle of the energy bins

    # Maxwellian spectra
    # π is gone as we do not normalize in /ster
    # However we keep the 2 as it seems necessary to keep the total energy flux
    # Either the Meier et al. conversion to ster is slightly wrong, or there is something I don't get?
    Φₘ = IeE_tot_eV / (2 * E₀^3) .* E_middle .* exp.(-E_middle ./ E₀)
    # Parameter for the LET (corrected equations to match the Fig. 4)
    b = (0.8 * E₀) .* (E₀ < 500) +
        (0.1 * E₀ + 350) .* (E₀ >= 500) # use 350 instead of .35 because we are in eV

    for i_μ in Beams
        if low_energy_tail
            Ie_max = maximum(Φₘ) # max of the Maxwellian - to scale LET amplitude
            Ie_top[i_μ, 1, :] = (Φₘ .+ 0.4 * Ie_max * (E₀ ./ E_middle) .* exp.(-E_middle ./ b)) .*
                                dE * BeamWeight[i_μ] ./ sum(BeamWeight[Beams])
        else
            Ie_top[i_μ, 1, :] = Φₘ .*
                                dE * BeamWeight[i_μ] ./ sum(BeamWeight[Beams])
        end
    end

    return Ie_top
end
