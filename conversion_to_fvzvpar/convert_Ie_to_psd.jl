using LoopVectorization

"""
    convert_Ie_to_psd.jl

Convert AURORA electron number flux `Ie` [#/m²/s] into electron phase-space density.

## Physics background

`Ie[j, i_z, i_t, i]` is the total flux (electrons per m² per second) **integrated** over
energy bin ΔE[i] and solid-angle bin ΔΩ[j].  The phase-space density is:

    f(v‖, v⊥) = Ie / (ΔE[i] * ΔΩ[j] * v[i]²)   [s³ m⁻⁶]

and the reduced (1-D marginal) distribution is:

    F(v‖) = ∫ f 2π v⊥ dv⊥  →  F[j,…,i] = Ie[j,…,i] / ΔE[i]   [s m⁻⁴]

Solid-angle bin widths (exact, no small-angle approximation):

    ΔΩ[j] = 2π * |cos(α_low[j]) - cos(α_high[j])|

where α is the pitch angle (0° = field-aligned up, 180° = field-aligned down).

## Coordinate convention

`v‖ = v * cos(α_center)`:
  - v‖ < 0  →  upward   (α < 90°)
  - v‖ > 0  →  downward (α > 90°)

`v⊥ = v * sin(α_center) ≥ 0`

## Units

| Quantity | Unit     |
|----------|----------|
| E        | eV       |
| v        | m/s      |
| Ie       | m⁻² s⁻¹ (integrated over ΔE and ΔΩ) |
| f        | s³ m⁻⁶   |
| F        | s m⁻⁴    |

## Exports

- `load_Ie(path)` — load a `IeFlickering-NN.mat` file and return named tuple
- `psd_grids(E_eV, mu_lims_cosine)` — compute velocity grids and bin widths
- `compute_f(Ie, ΔE_J, ΔΩ, v)` — full 2-D phase-space density [n_μ, n_z, n_t, n_E]
- `compute_F(Ie, ΔE_J)` — reduced 1-D distribution [n_μ, n_z, n_t, n_E]
- `make_psd_from_AURORA(path)` — top-level convenience wrapper
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
- `Ie`      : `[Nz, n_μ, Nt, nE]` — number flux [m⁻² s⁻¹], integrated over ΔE and ΔΩ
- `E`       : energy bin *left edges* [eV], length `nE`
- `mu_lims` : cosine pitch-angle bin *limits*, length `n_μ + 1`
- `t_run`   : time vector [s], length `Nt`
- `h_atm`   : altitude grid [m], length `Nz`
"""
function load_Ie(path::AbstractString)
    file = matopen(path)
        Ie_raw  = read(file, "Ie_ztE")   # [n_μ*Nz, Nt, nE]
        E       = vec(read(file, "E"))
        mu_lims = vec(read(file, "mu_lims"))
        t_run   = vec(read(file, "t_run"))
        h_atm   = vec(read(file, "h_atm"))
    close(file)

    nE  = size(Ie_raw, 3)
    Nz  = length(h_atm)
    n_μ = length(mu_lims) - 1
    Nt  = length(t_run)

    # Reshape from [n_μ*Nz, Nt, nE] → [Nz, n_μ, Nt, nE]
    Ie = reshape(Ie_raw, Nz, n_μ, Nt, nE)

    return (Ie=Ie, E=E[1:nE], mu_lims=mu_lims, t_run=t_run, h_atm=h_atm)
end

# ---------------------------------------------------------------------------
# Grid helpers
# ---------------------------------------------------------------------------

"""
    psd_grids(E_eV, mu_lims_cosine)

Compute velocity grids and bin widths from the AURORA energy and pitch-angle grids.

# Arguments
- `E_eV`           : energy left-edge grid [eV], length `nE`
- `mu_lims_cosine` : cosine of pitch-angle limits, length `n_μ + 1`

# Returns named tuple
- `v`         : speed at energy bin centres [m/s], length `nE`
- `ΔE_J`      : energy bin widths [J], length `nE`
- `ΔΩ`        : solid-angle bin widths [sr], length `n_μ`
- `μ_center`  : flux-weighted average cosine per pitch bin, length `n_μ`
                (= ⟨cosα⟩ with sinα weighting; see `mu_avg` in AURORA)
- `α_center`  : corresponding pitch-angle bin centres [rad], length `n_μ`
- `v_par`     : v‖ at every (pitch, energy) bin centre [m/s], matrix `[n_μ, nE]`
- `v_perp`    : v⊥ at every (pitch, energy) bin centre [m/s], matrix `[n_μ, nE]`
"""
function psd_grids(E_eV::AbstractVector, mu_lims_cosine::AbstractVector)
    nE  = length(E_eV)
    n_μ = length(mu_lims_cosine) - 1

    # --- Energy bin widths and bin-centre energies ---
    # E_eV holds the left edges of each bin; the right edge of bin i is E_eV[i+1]
    # (the last bin repeats the previous width)
    ΔE_eV = Vector{Float64}(undef, nE)
    for i in 1:(nE - 1)
        ΔE_eV[i] = E_eV[i + 1] - E_eV[i]
    end
    ΔE_eV[end] = ΔE_eV[end - 1]   # repeat last width

    E_centre_eV = E_eV .+ 0.5 .* ΔE_eV
    ΔE_J        = ΔE_eV .* e_charge

    # --- Speed at energy bin centres ---
    v = sqrt.(2 .* (E_centre_eV .* e_charge) ./ m_e)   # [m/s]

    # --- Solid-angle bin widths ---
    # Use abs() to be robust against either ascending or descending ordering of mu_lims.
    ΔΩ = Vector{Float64}(undef, n_μ)
    for j in 1:n_μ
        ΔΩ[j] = 2π * abs(mu_lims_cosine[j] - mu_lims_cosine[j + 1])
    end

    # --- Pitch-angle bin centres (flux-weighted) ---
    # mu_avg expects pitch-angle limits in degrees; returns ⟨cosα⟩ per bin.
    θ_lims_deg = acosd.(mu_lims_cosine)
    μ_center   = mu_avg(θ_lims_deg)                         # [n_μ]
    α_center   = acos.(clamp.(μ_center, -1.0, 1.0))         # [rad]

    # --- Velocity coordinates at each (j, i) bin centre ---
    v_par  = [μ_center[j] * v[i]        for j in 1:n_μ, i in 1:nE]   # [n_μ, nE]
    v_perp = [sin(α_center[j]) * v[i]   for j in 1:n_μ, i in 1:nE]

    # --- Parallel velocity-space width of each (j, i) bin ---
    # From dv = dE / (m_e * v), projected: dv‖ = |cos(α)| * dE / (m_e * v)
    dv_par = [abs(μ_center[j]) * ΔE_J[i] / (m_e * v[i]) for j in 1:n_μ, i in 1:nE]  # [n_μ, nE]

    return (v=v, ΔE_J=ΔE_J, ΔΩ=ΔΩ, μ_center=μ_center, α_center=α_center,
            v_par=v_par, v_perp=v_perp, dv_par=dv_par)
end

# ---------------------------------------------------------------------------
# Phase-space density — full 2-D  f(v‖, v⊥)
# ---------------------------------------------------------------------------

"""
    compute_f(Ie, ΔE_J, ΔΩ, v) -> f

Compute full phase-space density from AURORA electron flux.

# Arguments
- `Ie`   : `[Nz, n_μ, Nt, nE]` — flux [m⁻² s⁻¹], integrated over ΔE and ΔΩ
- `ΔE_J` : energy bin widths [J], length `nE`
- `ΔΩ`   : solid-angle bin widths [sr], length `n_μ`
- `v`    : speed at energy bin centres [m/s], length `nE`

# Returns
- `f`    : `[Nz, n_μ, Nt, nE]` — phase-space density [s³ m⁻⁶]

# Formula
    f[iz, j, it, i] = Ie[iz, j, it, i] / (ΔE_J[i] * ΔΩ[j] * v[i]²)
"""
function compute_f(Ie::AbstractArray{<:Real,4}, ΔE_J::AbstractVector,
                   ΔΩ::AbstractVector, v::AbstractVector)
    Nz, n_μ, Nt, nE = size(Ie)
    @assert length(ΔE_J) == nE  "ΔE_J length must match energy dimension of Ie"
    @assert length(ΔΩ)   == n_μ "ΔΩ length must match pitch-angle dimension of Ie"
    @assert length(v)    == nE  "v length must match energy dimension of Ie"

    f = similar(Ie, Float64)
    @tturbo for i in 1:nE
        denom_E = ΔE_J[i] * v[i]^2
        for j in 1:n_μ
            denom = denom_E * ΔΩ[j]
            for it in 1:Nt, iz in 1:Nz
                f[iz, j, it, i] = Ie[iz, j, it, i] / denom
            end
        end
    end
    return f
end

# ---------------------------------------------------------------------------
# Reduced 1-D distribution  F(v‖)  —  binned onto a regular v‖ grid
# ---------------------------------------------------------------------------

"""
    compute_F(Ie, ΔE_J, v_par, dv_par; n_vpar_bins=200, vpar_edges=nothing) -> F_binned, vpar_centres

Compute the reduced 1-D distribution F(v‖) binned onto a uniform v‖ grid,
conserving the integral ∫ F dv‖ and the total electron count.

# Arguments
- `Ie`         : `[Nz, n_μ, Nt, nE]` — flux [m⁻² s⁻¹], integrated over ΔE and ΔΩ
- `ΔE_J`       : energy bin widths [J], length `nE`
- `v_par`      : signed v‖ at each (j,i) bin centre [m/s], matrix `[n_μ, nE]`
- `dv_par`     : parallel velocity width of each (j,i) bin [m/s], matrix `[n_μ, nE]`
- `n_vpar_bins`: number of output v‖ bins (default 200); ignored if `vpar_edges` is given
- `vpar_edges` : custom bin edges [m/s], length `Nvpar + 1` (optional).
                 Must be a sorted AbstractVector. Input points outside this range are discarded.

# Returns
- `F_binned`     : `[Nvpar, Nz, Nt]` — F(v‖) [s m⁻⁴] on the output v‖ grid
- `vpar_centres` : `[Nvpar]` — v‖ bin centres [m/s]

# Physics
Each input point (j,i) represents a density F = Ie/ΔE [s m⁻⁴] spread uniformly over
the interval [v_par - dv‖/2, v_par + dv‖/2].  The fraction of that interval overlapping
each output bin b is computed exactly, so both ∫ F dv‖ and total electron count are
conserved regardless of whether the output edges align with the input bin edges.

Deposit into bin b:  F[j,…,i] * overlap(j,i,b)
Then divide by Δvbin[b] to recover density.  Non-uniform output grids are supported.
"""
function compute_F(Ie::AbstractArray{<:Real,4}, ΔE_J::AbstractVector,
                   v_par::AbstractMatrix, dv_par::AbstractMatrix;
                   n_vpar_bins::Int=200,
                   vpar_edges::Union{AbstractVector,Nothing}=nothing)
    Nz, n_μ, Nt, nE = size(Ie)
    @assert length(ΔE_J)  == nE        "ΔE_J must match energy dimension of Ie"
    @assert size(v_par)   == (n_μ, nE) "v_par must be [n_μ, nE]"
    @assert size(dv_par)  == (n_μ, nE) "dv_par must be [n_μ, nE]"

    # --- Build output grid ---
    edges = if isnothing(vpar_edges)
        range(minimum(v_par), maximum(v_par), length=n_vpar_bins + 1)
    else
        @assert issorted(vpar_edges) "vpar_edges must be sorted"
        vpar_edges
    end
    Nvpar        = length(edges) - 1
    vpar_centres = [(edges[b] + edges[b+1]) / 2 for b in 1:Nvpar]
    Δvbin        = [edges[b+1] - edges[b] for b in 1:Nvpar]
    inv_Δvbin    = 1.0 ./ Δvbin

    # --- Precompute overlap weights (independent of iz, it) ---
    # For each (j,i), find all output bins b that overlap with
    # [v_par[j,i] - dv‖/2, v_par[j,i] + dv‖/2] and store (b, overlap_length).
    # Deposit = F[j,…,i] * overlap_length;  then divide by Δvbin[b] at the end.
    overlaps = [Tuple{Int,Float64}[] for _ in 1:(n_μ * nE)]

    @inbounds for i in 1:nE, j in 1:n_μ
        dv   = dv_par[j, i]
        dv == 0.0 && continue
        v_lo = v_par[j, i] - 0.5 * dv
        v_hi = v_par[j, i] + 0.5 * dv
        # first candidate bin: last edge ≤ v_lo
        b_start = max(1, searchsortedlast(edges, v_lo))
        for b in b_start:Nvpar
            edges[b] >= v_hi && break
            ov = min(v_hi, edges[b+1]) - max(v_lo, edges[b])
            ov > 0.0 && push!(overlaps[(i-1)*n_μ + j], (b, ov))
        end
    end

    F_binned = zeros(Nvpar, Nz, Nt)

    Threads.@threads for it in 1:Nt
        for iz in 1:Nz
            @inbounds for i in 1:nE
                inv_dE = 1.0 / ΔE_J[i]
                for j in 1:n_μ
                    Ie_val = Ie[iz, j, it, i] * inv_dE   # F = Ie/ΔE [s m⁻⁴]
                    for (b, ov) in overlaps[(i-1)*n_μ + j]
                        F_binned[b, iz, it] += Ie_val * ov
                    end
                end
            end
            # divide each bin's accumulated integral by its width → density [s m⁻⁴]
            for b in 1:Nvpar
                F_binned[b, iz, it] *= inv_Δvbin[b]
            end
        end
    end

    return F_binned, vpar_centres
end

# ---------------------------------------------------------------------------
# Top-level convenience wrapper
# ---------------------------------------------------------------------------

"""
    make_psd_from_AURORA(path_to_file)

Load an AURORA result file and return both the full phase-space density `f(v‖, v⊥)` and
the reduced distribution `F(v‖)`, together with all grids needed for analysis and plotting.

# Arguments
- `path_to_file` : path to a `.mat` file containing `Ie_ztE`, `E`, `mu_lims`, `t_run`, `h_atm`
- `compute`      : what to compute — `:both` (default), `:f_only`, or `:F_only`

# Returns named tuple
- `f`            : `[Nz, n_μ, Nt, nE]` — full phase-space density [s³ m⁻⁶] (only if `compute ∈ (:both, :f_only)`)
- `F`            : `[Nvpar, Nz, Nt]`   — reduced distribution F(v‖) [s m⁻⁴] on uniform v‖ grid (only if `compute ∈ (:both, :F_only)`)
- `vpar_centres` : `[Nvpar]`           — v‖ bin centres [m/s]
- `v`            : speed grid [m/s], length `nE`
- `v_par`        : v‖ at each (pitch,energy) bin [m/s], matrix `[n_μ, nE]`
- `v_perp`       : v⊥ at each (pitch,energy) bin [m/s], matrix `[n_μ, nE]`
- `dv_par`       : Δv‖ width of each (pitch,energy) bin [m/s], matrix `[n_μ, nE]`
- `ΔE_J`         : energy bin widths [J], length `nE`
- `ΔΩ`           : solid-angle bin widths [sr], length `n_μ`
- `μ_center`     : average cosine per pitch bin, length `n_μ`
- `E`            : energy left-edge grid [eV], length `nE`
- `t_run`        : time vector [s]
- `h_atm`        : altitude grid [m]
- `mu_lims`      : cosine pitch-angle limits, length `n_μ + 1`

# Example
```julia
include("convert_Ie_to_psd.jl")
res = make_psd_from_AURORA("data/IeFlickering-01.mat")              # both f and F
res = make_psd_from_AURORA("data/IeFlickering-01.mat"; compute=:f_only)  # only f
res = make_psd_from_AURORA("data/IeFlickering-01.mat"; compute=:F_only)  # only F
# res.f[iz, j, it, i]  → phase-space density at alt iz, beam j, time it, energy i
# res.v_par[j, i]       → corresponding v‖ coordinate [m/s]
```
"""
function make_psd_from_AURORA(path_to_file::AbstractString;
                               compute::Symbol=:both,
                               n_vpar_bins::Int=200,
                               vpar_edges::Union{AbstractVector,Nothing}=nothing)
    compute ∈ (:both, :f_only, :F_only) ||
        throw(ArgumentError("compute must be :both, :f_only, or :F_only; got :$compute"))

    data  = load_Ie(path_to_file)
    grids = psd_grids(data.E, data.mu_lims)

    f = if compute ∈ (:both, :f_only)
        compute_f(data.Ie, grids.ΔE_J, grids.ΔΩ, grids.v)
    else
        nothing
    end

    F_binned, vpar_centres = if compute ∈ (:both, :F_only)
        compute_F(data.Ie, grids.ΔE_J, grids.v_par, grids.dv_par;
                  n_vpar_bins, vpar_edges)
    else
        nothing, nothing
    end

    return (f            = f,
            F            = F_binned,       # [Nvpar, Nz, Nt]  F(v‖) [s m⁻⁴]
            vpar_centres = vpar_centres,   # [Nvpar]           v‖ bin centres [m/s]
            v            = grids.v,
            v_par        = grids.v_par,
            v_perp       = grids.v_perp,
            dv_par       = grids.dv_par,
            ΔE_J         = grids.ΔE_J,
            ΔΩ           = grids.ΔΩ,
            μ_center     = grids.μ_center,
            E            = data.E,
            t_run        = data.t_run,
            h_atm        = data.h_atm,
            mu_lims      = data.mu_lims)
end
