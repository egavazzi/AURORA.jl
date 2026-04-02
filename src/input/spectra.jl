## ====================================================================================== ##
## Abstract type
## ====================================================================================== ##

"""
    AbstractSpectrum

Abstract supertype for all electron energy spectrum shapes.

Concrete subtypes: [`FlatSpectrum`](@ref), [`GaussianSpectrum`](@ref),
[`MaxwellianSpectrum`](@ref), [`FileSpectrum`](@ref).
"""
abstract type AbstractSpectrum end



## ====================================================================================== ##
## FlatSpectrum
## ====================================================================================== ##

"""
    FlatSpectrum(IeE_tot; E_min=0.0)

Flat (constant) differential number flux spectrum above a minimum energy.

The flux is uniform in #e⁻/m²/s/eV for energies above `E_min`, and zero below.
Normalized so that the integrated energy flux equals `IeE_tot`.

# Arguments
- `IeE_tot`: total energy flux (W/m²)

# Keyword Arguments
- `E_min=0.0`: minimum energy threshold (eV). Flux is zero below this energy.

# Examples
```julia
FlatSpectrum(1e-2)                  # flat from lowest energy
FlatSpectrum(1e-2; E_min=2900.0)    # flat above 2900 eV
```
"""
struct FlatSpectrum <: AbstractSpectrum
    IeE_tot::Float64
    E_min::Float64

    function FlatSpectrum(IeE_tot; E_min=0.0)
        IeE_tot > 0 || error("IeE_tot must be positive, got $IeE_tot")
        E_min >= 0 || error("E_min must be non-negative, got $E_min")
        new(Float64(IeE_tot), Float64(E_min))
    end
end

function Base.show(io::IO, s::FlatSpectrum)
    print(io, "FlatSpectrum(IeE_tot=$(s.IeE_tot) W/m², E_min=$(s.E_min) eV)")
end



## ====================================================================================== ##
## GaussianSpectrum
## ====================================================================================== ##

"""
    GaussianSpectrum(IeE_tot, E₀, ΔE)

Gaussian energy spectrum centered at `E₀` with width `ΔE`.

The spectral shape is Φ(E) ∝ exp(-(E - E₀)² / ΔE²), normalized so that the
integrated energy flux equals `IeE_tot`.

# Arguments
- `IeE_tot`: total energy flux (W/m²)
- `E₀`: center energy (eV)
- `ΔE`: energy width (eV)

# Examples
```julia
GaussianSpectrum(1e-2, 5000.0, 500.0)
```
"""
struct GaussianSpectrum <: AbstractSpectrum
    IeE_tot::Float64
    E₀::Float64
    ΔE::Float64

    function GaussianSpectrum(IeE_tot, E₀, ΔE)
        IeE_tot > 0 || error("IeE_tot must be positive, got $IeE_tot")
        E₀ > 0 || error("E₀ must be positive, got $E₀")
        ΔE > 0 || error("ΔE must be positive, got $ΔE")
        new(Float64(IeE_tot), Float64(E₀), Float64(ΔE))
    end
end

function Base.show(io::IO, s::GaussianSpectrum)
    print(io, "GaussianSpectrum(IeE_tot=$(s.IeE_tot) W/m², E₀=$(s.E₀) eV, ΔE=$(s.ΔE) eV)")
end



## ====================================================================================== ##
## MaxwellianSpectrum
## ====================================================================================== ##

"""
    MaxwellianSpectrum(IeE_tot, E₀; low_energy_tail=true)

Maxwellian energy spectrum with optional low-energy tail (LET).

Based on the corrected implementation of Meier/Strickland/Hecht/Christensen
JGR 1989 (pages 13541-13552).

# Arguments
- `IeE_tot`: total energy flux (W/m²)
- `E₀`: characteristic energy (eV)

# Keyword Arguments
- `low_energy_tail=true`: include a low-energy tail following Meier et al. 1989

# Examples
```julia
MaxwellianSpectrum(1e-2, 1000.0)
MaxwellianSpectrum(1e-2, 1000.0; low_energy_tail=false)
```
"""
struct MaxwellianSpectrum <: AbstractSpectrum
    IeE_tot::Float64
    E₀::Float64
    low_energy_tail::Bool

    function MaxwellianSpectrum(IeE_tot, E₀; low_energy_tail=true)
        IeE_tot > 0 || error("IeE_tot must be positive, got $IeE_tot")
        E₀ > 0 || error("E₀ must be positive, got $E₀")
        new(Float64(IeE_tot), Float64(E₀), low_energy_tail)
    end
end

function Base.show(io::IO, s::MaxwellianSpectrum)
    tail = s.low_energy_tail ? " + LET" : ""
    print(io, "MaxwellianSpectrum(IeE_tot=$(s.IeE_tot) W/m², E₀=$(s.E₀) eV$(tail))")
end



## ====================================================================================== ##
## FileSpectrum
## ====================================================================================== ##

"""
    FileSpectrum(filename; interpolation=:constant)

Energy spectrum loaded from a `.mat` file.

The file contains the complete flux distribution (possibly time-dependent).
When used inside an [`InputFlux`](@ref), this spectrum type requires
[`ConstantModulation`](@ref) because the file already contains the full temporal behavior.

# Arguments
- `filename`: path to the `.mat` file containing the flux data

# Keyword Arguments
- `interpolation=:constant`: interpolation scheme for resampling the file's time grid.
  Either `:constant`, `:linear`, or `:pchip`.

# File format
The `.mat` file must contain:
- `Ie_total`: flux array of shape `[n_μ, n_t_file, n_E]`
- `t_top`: time grid (s), as a vector `[n_t_file]`. Optional when `n_t_file == 1`.

# Examples
```julia
FileSpectrum("path/to/flux.mat")
FileSpectrum("path/to/flux.mat"; interpolation=:linear)
```
"""
struct FileSpectrum <: AbstractSpectrum
    filename::String
    interpolation::Symbol

    function FileSpectrum(filename; interpolation=:constant)
        interpolation in (:constant, :linear, :pchip) ||
            error("Unknown interpolation type: $interpolation. Must be :constant, :linear, or :pchip.")
        new(filename, interpolation)
    end
end

function Base.show(io::IO, s::FileSpectrum)
    print(io, "FileSpectrum(\"$(basename(s.filename))\", interpolation=$(s.interpolation))")
end



## ====================================================================================== ##
## evaluate_spectrum
## ====================================================================================== ##

"""
    evaluate_spectrum(spec::AbstractSpectrum, model::AuroraModel)

Evaluate the energy spectrum on the model's energy grid, returning a vector of
differential number flux per eV (#e⁻/m²/s/eV) of length `n_E`.

This function does **not** distribute over beams or apply temporal modulation.
Use [`compute_flux`](@ref) for the full 3D flux array.
"""
function evaluate_spectrum end


"""
    evaluate_spectrum(spec::FlatSpectrum, model::AuroraModel)

Evaluate a flat spectrum on the model's energy grid.
"""
function evaluate_spectrum(spec::FlatSpectrum, model::AuroraModel)
    E_centers = model.energy_grid.E_centers
    ΔE = model.energy_grid.ΔE

    if spec.E_min > E_centers[end]
        error("E_min ($(spec.E_min) eV) is larger than the maximum energy in the grid ($(E_centers[end]) eV).")
    end

    qₑ = 1.602176620898e-19
    IeE_tot_eV = spec.IeE_tot / qₑ
    return _flat_spectrum(IeE_tot_eV, E_centers, ΔE, spec.E_min)
end


"""
    evaluate_spectrum(spec::GaussianSpectrum, model::AuroraModel)

Evaluate a Gaussian spectrum on the model's energy grid.
"""
function evaluate_spectrum(spec::GaussianSpectrum, model::AuroraModel)
    E_centers = model.energy_grid.E_centers
    ΔE = model.energy_grid.ΔE

    qₑ = 1.602176620898e-19
    IeE_tot_eV = spec.IeE_tot / qₑ
    return _gaussian_spectrum(IeE_tot_eV, E_centers, ΔE, spec.E₀, spec.ΔE)
end


"""
    evaluate_spectrum(spec::MaxwellianSpectrum, model::AuroraModel)

Evaluate a Maxwellian spectrum (with optional low-energy tail) on the model's energy grid.
"""
function evaluate_spectrum(spec::MaxwellianSpectrum, model::AuroraModel)
    E_centers = model.energy_grid.E_centers
    ΔE = model.energy_grid.ΔE
    E₀ = spec.E₀

    qₑ = 1.602176620898e-19
    IeE_tot_eV = spec.IeE_tot / qₑ

    # Maxwellian spectra (corrected Meier et al. 1989)
    Φₘ = IeE_tot_eV / (2 * E₀^3) .* E_centers .* exp.(-E_centers ./ E₀)

    if spec.low_energy_tail
        # Parameter for the LET (corrected equations to match Fig. 4)
        b = (0.8 * E₀) * (E₀ < 500) + (0.1 * E₀ + 350) * (E₀ >= 500)
        Ie_max = maximum(Φₘ)
        Φ = Φₘ .+ 0.4 * Ie_max * (E₀ ./ E_centers) .* exp.(-E_centers ./ b)
    else
        Φ = Φₘ
    end

    return Φ  # #e⁻/m²/s/eV (differential, before multiplying by ΔE)
end
