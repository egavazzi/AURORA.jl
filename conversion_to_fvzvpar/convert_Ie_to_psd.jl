using LoopVectorization

"""
    convert_Ie_to_psd.jl

Convert AURORA electron number flux `Ie` [#/m²/s] into electron phase-space density
`f(v_parallel, v_perp)` [s³ m⁻⁶] on the original `(E, Ω)` bin layout.

The output from AURORA `Ie[iz, j, it, i]` is already integrated over energy bin
`ΔE[i]` and solid-angle bin `ΔΩ[j]`, so the conversion is:

    f = Ie * m_e / (ΔE * ΔΩ * v^2)

where `v` is the speed at the energy-bin centre.
"""

using MAT: matopen
using AURORA: mu_avg

const m_e      = 9.10938356e-31   # kg
const e_charge = 1.602176634e-19  # C

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

"""
    load_Ie(path)

Load a `IeFlickering-NN.mat` result file produced by AURORA.

Returns a named tuple with fields:
- `Ie`      : `[Nz, n_mu, Nt, nE]` number flux [m⁻² s⁻¹], integrated over ΔE and ΔΩ
- `E`       : energy bin left edges [eV], length `nE`
- `mu_lims` : cosine pitch-angle bin limits, length `n_mu + 1`
- `t_run`   : time vector [s], length `Nt`
- `h_atm`   : altitude grid [m], length `Nz`
"""
function load_Ie(path::AbstractString)
    file = matopen(path)
    Ie_raw  = read(file, "Ie_ztE")   # [n_mu*Nz, Nt, nE]
    E       = vec(read(file, "E"))
    μ_lims = vec(read(file, "mu_lims"))
    t_run   = vec(read(file, "t_run"))
    h_atm   = vec(read(file, "h_atm"))
    close(file)

    nE  = size(Ie_raw, 3)
    Nz  = length(h_atm)
    nμ = length(μ_lims) - 1
    Nt  = length(t_run)

    # Reshape from [n_mu*Nz, Nt, nE] -> [Nz, n_mu, Nt, nE]
    Ie = reshape(Ie_raw, Nz, nμ, Nt, nE)

    return (Ie=Ie, E=E[1:nE], μ_lims=μ_lims, t_run=t_run, h_atm=h_atm)
end

# ---------------------------------------------------------------------------
# Grid helpers for Ie -> f
# ---------------------------------------------------------------------------

"""
    psd_grids(E_eV, μ_lims_cosine)

Compute energy, angle, and velocity grids needed by `compute_f`.

Returns a named tuple with:
- `v`        : speed at energy-bin centres [m/s], length `nE`
- `ΔE_J` : energy-bin widths [J], length `nE`
- `ΔΩ` : solid-angle widths [sr], length `nμ`
- `μ_center`  : flux-weighted average cosine per pitch bin, length `nμ`
- `α_center` : pitch-angle centres [rad], length `nμ`
- `v_par`    : `v_parallel` at each (pitch, energy) bin centre [m/s], `[nμ, nE]`
- `v_perp`   : `v_perp` at each (pitch, energy) bin centre [m/s], `[nμ, nE]`
"""
function psd_grids(E_eV::AbstractVector, μ_lims_cosine::AbstractVector)
    nE   = length(E_eV)
    nμ = length(μ_lims_cosine) - 1

    ΔE_eV = diff(E_eV); ΔE_eV = vcat(ΔE_eV, ΔE_eV[end])
    ΔE_J    = ΔE_eV .* e_charge
    E_center_eV = E_eV .+ 0.5 .* ΔE_eV
    E_center_J = E_center_eV .* e_charge
    v           = sqrt.(2 .* E_center_J ./ m_e)

    ΔΩ = Vector{Float64}(undef, nμ)
    for j in 1:nμ
        ΔΩ[j] = 2pi * abs(μ_lims_cosine[j] - μ_lims_cosine[j + 1])
    end

    theta_lims_deg = acosd.(μ_lims_cosine)
    μ_center      = mu_avg(theta_lims_deg)
    α_center   = acos.(clamp.(μ_center, -1.0, 1.0))

    v_par  = [μ_center[j] * v[i]      for j in 1:nμ, i in 1:nE]
    v_perp = [sin(α_center[j]) * v[i] for j in 1:nμ, i in 1:nE]

    return (v=v,
            ΔE_J=ΔE_J,
            ΔΩ=ΔΩ,
            μ_center=μ_center,
            α_center=α_center,
            v_par=v_par,
            v_perp=v_perp)
end

# ---------------------------------------------------------------------------
# Ie -> f conversion
# ---------------------------------------------------------------------------

"""
    compute_f(Ie, ΔE_J, ΔΩ, v) -> f

Convert AURORA electron flux `Ie` to full phase-space density `f` on the same
`[Nz, nμ, Nt, nE]` indexing.

Formula:

    f[iz, j, it, i] = Ie[iz, j, it, i] / (ΔE_J[i] * ΔΩ[j]) * m_e / v[i]^2
"""
function compute_f(Ie::AbstractArray{<:Real,4}, ΔE_J::AbstractVector,
                   ΔΩ::AbstractVector, v::AbstractVector)
    Nz, nμ, Nt, nE = size(Ie)
    @assert length(ΔE_J) == nE  "ΔE_J length must match energy dimension of Ie"
    @assert length(ΔΩ) == nμ "ΔΩ length must match pitch-angle dimension of Ie"
    @assert length(v) == nE "v length must match energy dimension of Ie"

    f = similar(Ie, Float64)
    @tturbo for i in 1:nE
        denom_E = ΔE_J[i] * v[i]^2
        for j in 1:nμ
            denom = denom_E * ΔΩ[j]
            for it in 1:Nt, iz in 1:Nz
                f[iz, j, it, i] = Ie[iz, j, it, i] * m_e / denom
            end
        end
    end
    return f
end

# ---------------------------------------------------------------------------
# Ie -> F(v_parallel) reduction
# ---------------------------------------------------------------------------

"""
    default_vpar_edges(v, μ_lims) -> vpar_edges

Construct a signed `v_parallel` grid from all source-bin endpoints `μ_lims * v`.
This guarantees that the reduction grid spans the populated phase-space region
without truncating the outermost bins.
"""
function default_vpar_edges(v::AbstractVector, μ_lims::AbstractVector)
    vpar_edges = sort!(vec([μ * vi for μ in μ_lims, vi in v]))
    unique!(vpar_edges)
    return vpar_edges
end

"""
    compute_F(Ie, μ_lims, v; vpar_edges=nothing) -> (F, vpar_edges, vpar_centres, Δvpar)

Reduce AURORA electron flux `Ie` to a 1D distribution `F(v_parallel)` while
conserving total electron density:

    sum(Ie / v) == integral(F dv_parallel)

Each original `(μ, E)` bin contributes the density packet `Ie / v` distributed
uniformly over the source interval `[μ_lo * v, μ_hi * v]` and deposited onto the
target `v_parallel` bins by interval overlap.
"""
function compute_F(Ie::AbstractArray{<:Real,4}, μ_lims::AbstractVector,
                   v::AbstractVector; vpar_edges::Union{Nothing,AbstractVector}=nothing)
    Nz, nμ, Nt, nE = size(Ie)
    @assert length(μ_lims) == nμ + 1 "μ_lims length must match pitch-angle dimension of Ie"
    @assert length(v) == nE "v length must match energy dimension of Ie"

    if isnothing(vpar_edges)
        vpar_edges = default_vpar_edges(v, μ_lims)
    else
        vpar_edges = collect(Float64, vpar_edges)
    end

    @assert length(vpar_edges) >= 2 "vpar_edges must contain at least two entries"
    @assert issorted(vpar_edges) "vpar_edges must be sorted in ascending order"

    Δvpar = diff(vpar_edges)
    @assert all(Δvpar .> 0) "vpar_edges must be strictly increasing"

    Nvpar = length(Δvpar)
    vpar_centres = @. 0.5 * (vpar_edges[1:end-1] + vpar_edges[2:end])
    F = zeros(Float64, Nvpar, Nz, Nt)

    for i in 1:nE
        inv_v = 1 / v[i]
        for j in 1:nμ
            src_lo = μ_lims[j] * v[i]
            src_hi = μ_lims[j + 1] * v[i]
            if src_lo > src_hi
                src_lo, src_hi = src_hi, src_lo
            end

            src_width = src_hi - src_lo
            if iszero(src_width)
                continue
            end

            k_start = max(1, searchsortedlast(vpar_edges, src_lo))
            k_end = min(Nvpar, searchsortedfirst(vpar_edges, src_hi) - 1)

            for k in k_start:k_end
                overlap_lo = max(src_lo, vpar_edges[k])
                overlap_hi = min(src_hi, vpar_edges[k + 1])
                overlap = overlap_hi - overlap_lo
                if overlap <= 0
                    continue
                end

                weight = overlap / src_width / Δvpar[k]
                @inbounds for it in 1:Nt, iz in 1:Nz
                    F[k, iz, it] += Ie[iz, j, it, i] * inv_v * weight
                end
            end
        end
    end

    return (F=F,
            vpar_edges=vpar_edges,
            vpar_centres=vpar_centres,
            Δvpar=Δvpar)
end

# ---------------------------------------------------------------------------
# Top-level convenience wrapper
# ---------------------------------------------------------------------------

"""
    make_psd_from_AURORA(path_to_file; compute=:f_only, vpar_edges=nothing)

Load an AURORA result file and compute full phase-space density `f`, the reduced
distribution `F(v_parallel)`, or both.

Returns a named tuple with:
- `f`            : `[Nz, nμ, Nt, nE]` phase-space density [s^3 m^-6]
- `F`            : `[Nvpar, Nz, Nt]` reduced distribution [s m^-4]
- `vpar_edges`   : `v_parallel` bin edges [m/s], length `Nvpar + 1`
- `vpar_centres` : `v_parallel` bin centres [m/s], length `Nvpar`
- `Δvpar`        : `v_parallel` bin widths [m/s], length `Nvpar`
- `v`            : speed grid [m/s], length `nE`
- `v_par`        : `v_parallel` at each (pitch,energy) bin centre [m/s], `[nμ, nE]`
- `v_perp`       : `v_perp` at each (pitch,energy) bin centre [m/s], `[nμ, nE]`
- `ΔE_J`     : energy-bin widths [J], length `nE`
- `ΔΩ`   : solid-angle widths [sr], length `nμ`
- `μ_center`    : average cosine per pitch bin, length `nμ`
- `E`            : energy left-edge grid [eV], length `nE`
- `t_run`        : time vector [s]
- `h_atm`        : altitude grid [m]
- `μ_lims`      : cosine pitch-angle limits, length `nμ + 1`
"""
function make_psd_from_AURORA(path_to_file::AbstractString;
                              compute::Symbol=:f_only,
                              vpar_edges::Union{Nothing,AbstractVector}=nothing)
    data  = load_Ie(path_to_file)
    grids = psd_grids(data.E, data.μ_lims)

    if compute ∉ (:f_only, :F_only, :both)
        throw(ArgumentError("compute must be one of :f_only, :F_only, or :both"))
    end

    f_result = compute in (:f_only, :both) ?
        (f=compute_f(data.Ie, grids.ΔE_J, grids.ΔΩ, grids.v),) :
        NamedTuple()

    F_result = compute in (:F_only, :both) ?
        compute_F(data.Ie, data.μ_lims, grids.v; vpar_edges=vpar_edges) :
        NamedTuple()

    return merge(f_result,
                 F_result,
                 (v          = grids.v,
                  v_par      = grids.v_par,
                  v_perp     = grids.v_perp,
                  ΔE_J       = grids.ΔE_J,
                  ΔΩ         = grids.ΔΩ,
                  μ_center   = grids.μ_center,
                  E          = data.E,
                  t_run      = data.t_run,
                  h_atm      = data.h_atm,
                  μ_lims     = data.μ_lims))
end
