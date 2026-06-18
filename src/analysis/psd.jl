using NCDatasets: NCDataset, defDim, defVar
using LoopVectorization
using ProgressMeter: Progress, next!

const m_e = 9.10938356e-31
const e_charge = 1.602176634e-19

"""
    compute_f(Ie, ΔE_J, ΔΩ, v) -> f

Convert AURORA electron flux `Ie` to full phase-space density `f` on the same
`[Nz, nμ, Nt, nE]` indexing.
"""
function compute_f(
    Ie::AbstractArray{<:Real, 4},
    ΔE_J::AbstractVector,
    ΔΩ::AbstractVector,
    v::AbstractVector,
)
    Nz, nμ, Nt, nE = size(Ie)
    @assert length(ΔE_J) == nE "ΔE_J length must match energy dimension of Ie"
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

"""
    compute_F(Ie, μ_lims, v; vpar_edges=nothing) -> (F, vpar_edges, vpar_centers, Δvpar)

Reduce AURORA electron flux `Ie` to a 1D distribution `F(v_parallel)` while
conserving total electron density.

If `vpar_edges` is `nothing`, a symmetric uniform grid is generated automatically
with `default_vpar_edges(v)`, spanning `[-maximum(v), maximum(v)]` and including
an exact edge at `v_parallel = 0`.
"""
function compute_F(
    Ie::AbstractArray{<:Real, 4},
    μ_lims::AbstractVector,
    v::AbstractVector;
    vpar_edges::Union{Nothing, AbstractVector} = nothing,
)
    Nz, nμ, Nt, nE = size(Ie)
    @assert length(μ_lims) == nμ + 1 "μ_lims length must match pitch-angle dimension of Ie"
    @assert length(v) == nE "v length must match energy dimension of Ie"

    if isnothing(vpar_edges)
        vpar_edges = default_vpar_edges(v)
    else
        vpar_edges = collect(Float64, vpar_edges)
    end

    @assert length(vpar_edges) >= 2 "vpar_edges must contain at least two entries"
    @assert issorted(vpar_edges) "vpar_edges must be sorted in ascending order"

    Δvpar = diff(vpar_edges)
    @assert all(Δvpar .> 0) "vpar_edges must be strictly increasing"

    Nvpar = length(Δvpar)
    vpar_centers = @. 0.5 * (vpar_edges[1:(end - 1)] + vpar_edges[2:end])
    F_local = zeros(Float64, Nz, Nvpar, Nt)

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
                @tturbo for it in 1:Nt, iz in 1:Nz
                    F_local[iz, k, it] += Ie[iz, j, it, i] * inv_v * weight
                end
            end
        end
    end

    F = permutedims(F_local, (2, 1, 3))

    return (F = F, vpar_edges = vpar_edges, vpar_centers = vpar_centers, Δvpar = Δvpar)
end


# ======================================================================================== #
#                                  I/O AND WRAPPERS                                      #
# ======================================================================================== #

"""
    make_psd_file(directory_to_process; compute=:both, vpar_edges=nothing, max_bytes=512*1024^2, compress=false, show_progress=false)

Read `simulation_data.nc` from `directory_to_process`, convert electron flux to
phase-space density, and write results to `analysis/psd.nc`.

The flux is streamed over time-chunks and the phase-space density is written to disk chunk
by chunk, so peak memory stays bounded even for large runs.

# Keyword Arguments
- `compute`: one of `:f_only`, `:F_only`, or `:both`.
- `vpar_edges`: custom `v_parallel` bin edges [m/s]. If `nothing`, an automatic
  symmetric uniform grid is used spanning `[-maximum(v), maximum(v)]`.
- `max_bytes`: per-chunk memory budget for streaming the flux (default 512 MiB).
- `compress`: zlib compression level for the `f` and `F` variables, with the same semantics
  as in [`AuroraOutputManager`](@ref): `false`/`0` (default, no compression), `true`
  (level 4), or an exact level `1`–`9`.
- `show_progress`: show a `ProgressMeter` progress bar while streaming chunks (default `false`).
"""
function make_psd_file(
    directory_to_process;
    compute::Symbol = :both,
    vpar_edges::Union{Nothing, AbstractVector} = nothing,
    max_bytes::Real = 512 * 1024^2,
    compress = false,
    show_progress::Bool = false,
)
    if compute ∉ (:f_only, :F_only, :both)
        throw(ArgumentError("compute must be one of :f_only, :F_only, or :both"))
    end
    dl = resolve_deflatelevel(compress)
    println("Converting Ie to PSD.")

    coord  = load_coordinates(directory_to_process)
    grids = psd_grids(coord.E_centers, coord.ΔE, coord.μ_lims)

    Nz, nμ, Nt, nE = coord.n_z, coord.n_μ, coord.n_t, coord.n_E
    want_f = compute in (:f_only, :both)
    want_F = compute in (:F_only, :both)

    # Fix the v_parallel grid up front so every time-chunk maps onto the same bins.
    vpe = want_F ?
          (isnothing(vpar_edges) ? default_vpar_edges(grids.v) : collect(Float64, vpar_edges)) :
          nothing

    analysis_dir = joinpath(directory_to_process, "analysis")
    mkpath(analysis_dir)
    savefile = joinpath(analysis_dir, "psd.nc")

    NCDataset(savefile, "c") do ds
        # Dimensions
        defDim(ds, "altitude",           Nz)
        defDim(ds, "time",               Nt)
        defDim(ds, "energy",             nE)
        defDim(ds, "pitch_angle",        nμ)
        defDim(ds, "pitch_angle_bounds", nμ + 1)
        defDim(ds, "energy_bounds",      nE + 1)

        # Coordinate variables
        alt_v = defVar(ds, "altitude", Float64, ("altitude",);
                       attrib=["units" => "m", "long_name" => "altitude"])
        alt_v[:] = coord.h_atm
        t_v = defVar(ds, "time", Float64, ("time",);
                     attrib=["units" => "s", "long_name" => "simulation time"])
        t_v[:] = coord.t
        en_v = defVar(ds, "energy", Float64, ("energy",);
                      attrib=["units" => "eV", "long_name" => "energy bin center"])
        en_v[:] = coord.E_centers
        ee_v = defVar(ds, "energy_edges", Float64, ("energy_bounds",);
                      attrib=["units" => "eV", "long_name" => "energy bin edges"])
        ee_v[:] = coord.E_edges
        pa_v = defVar(ds, "pitch_angle", Float64, ("pitch_angle",);
                      attrib=["units" => "1", "long_name" => "cosine of pitch angle (beam center)"])
        pa_v[:] = grids.μ_center
        ml_v = defVar(ds, "mu_lims", Float64, ("pitch_angle_bounds",);
                      attrib=["units" => "1", "long_name" => "pitch-angle cosine bin boundaries"])
        ml_v[:] = coord.μ_lims
        v_v = defVar(ds, "v", Float64, ("energy",);
                     attrib=["units" => "m s-1", "long_name" => "electron speed"])
        v_v[:] = grids.v
        bw_v = defVar(ds, "beam_weight", Float64, ("pitch_angle",);
                      attrib=["units" => "sr", "long_name" => "solid-angle beam weight"])
        bw_v[:] = grids.ΔΩ
        dE_v = defVar(ds, "dE_J", Float64, ("energy",);
                      attrib=["units" => "J", "long_name" => "energy bin width in Joules"])
        dE_v[:] = grids.ΔE_J

        # Data variables, filled by streaming below
        f_v = nothing
        if want_f
            f_v = defVar(ds, "f", Float64, ("altitude", "pitch_angle", "time", "energy");
                         deflatelevel=dl, shuffle=(dl > 0),
                         chunksizes=(Nz, nμ, 1, nE),
                         attrib=["units" => "s3 m-6", "long_name" => "phase-space density"])
        end

        F_v = nothing
        if want_F
            Nvpar = length(vpe) - 1
            defDim(ds, "vpar",        Nvpar)
            defDim(ds, "vpar_bounds", Nvpar + 1)
            vpc_v = defVar(ds, "vpar_centers", Float64, ("vpar",);
                           attrib=["units" => "m s-1", "long_name" => "v_parallel bin centers"])
            vpc_v[:] = @. 0.5 * (vpe[1:(end - 1)] + vpe[2:end])
            vpe_v = defVar(ds, "vpar_edges", Float64, ("vpar_bounds",);
                           attrib=["units" => "m s-1", "long_name" => "v_parallel bin edges"])
            vpe_v[:] = vpe
            dv_v = defVar(ds, "dvpar", Float64, ("vpar",);
                          attrib=["units" => "m s-1", "long_name" => "v_parallel bin widths"])
            dv_v[:] = diff(vpe)
            F_v = defVar(ds, "F", Float64, ("vpar", "altitude", "time");
                         deflatelevel=dl, shuffle=(dl > 0),
                         attrib=["units" => "s m-4",
                                  "long_name" => "reduced distribution function F(v_parallel)"])
        end

        # Stream the flux and fill the data variables chunk by chunk
        slice_bytes = Nz * nμ * nE * sizeof(Float64)
        nt_chunk = time_chunk_length(slice_bytes, max_bytes, Nt)
        n_chunks = ceil(Int, Nt / nt_chunk)
        progress = show_progress ? Progress(n_chunks; desc="Computing PSD: ", dt=1.0) : nothing
        foreach_Ie_time_chunk(directory_to_process; max_bytes) do Ie_chunk, t_range
            if want_f
                f_v[:, :, t_range, :] = compute_f(Ie_chunk, grids.ΔE_J, grids.ΔΩ, grids.v)
            end
            if want_F
                F_v[:, :, t_range] = compute_F(Ie_chunk, coord.μ_lims, grids.v; vpar_edges = vpe).F
            end
            isnothing(progress) || next!(progress)
        end
    end

    println("PSD saved in $savefile")
    return nothing
end

"""
    make_psd_file(sim::AuroraSimulation; kwargs...)

Convenience wrapper that calls [`make_psd_file`](@ref) on `sim.output.savedir`.
"""
make_psd_file(sim::AuroraSimulation; kwargs...) = make_psd_file(sim.output.savedir; kwargs...)


"""
    psd_grids(E_centers, ΔE, μ_lims_cosine)

Compute velocity and pitch-angle grids needed by `compute_f` and `compute_F`.

- `E_centers` is the energy bin centers (eV). Vector [nE]
- `ΔE` is the energy bin widths (eV). Vector [nE]
"""
function psd_grids(E_centers::AbstractVector, ΔE::AbstractVector, μ_lims_cosine::AbstractVector)
    nE = length(E_centers)
    nμ = length(μ_lims_cosine) - 1
    @assert length(ΔE) == nE "ΔE length must match E_centers length"

    ΔE_J = ΔE .* e_charge
    E_center_J = E_centers .* e_charge
    v = sqrt.(2 .* E_center_J ./ m_e)

    ΔΩ = Vector{Float64}(undef, nμ)
    for j in 1:nμ
        ΔΩ[j] = 2pi * abs(μ_lims_cosine[j] - μ_lims_cosine[j + 1])
    end

    theta_lims_deg = acosd.(μ_lims_cosine)
    μ_center = mu_avg(theta_lims_deg)
    α_center = acos.(clamp.(μ_center, -1.0, 1.0))

    v_par = [μ_center[j] * v[i] for j in 1:nμ, i in 1:nE]
    v_perp = [sin(α_center[j]) * v[i] for j in 1:nμ, i in 1:nE]

    return (
        v = v,
        ΔE_J = ΔE_J,
        ΔΩ = ΔΩ,
        μ_center = μ_center,
        α_center = α_center,
        v_par = v_par,
        v_perp = v_perp,
    )
end

"""
    default_vpar_edges(v) -> vpar_edges

Construct a symmetric, uniformly spaced `v_parallel` grid spanning
`[-maximum(v), maximum(v)]` with an exact edge at `v_parallel = 0`.
"""
function default_vpar_edges(v::AbstractVector; n_edges::Int = 75)
    vmax = maximum(v)
    positive_edges = collect(range(0.0, vmax; length = n_edges))
    return vcat(-reverse(positive_edges[2:end]), positive_edges)
end
