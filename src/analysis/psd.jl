using NCDatasets: NCDataset, defDim, defVar
using LoopVectorization

const m_e = 9.10938356e-31
const e_charge = 1.602176634e-19

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
    default_vpar_edges(v) -> vpar_edges

Construct a symmetric, uniformly spaced `v_parallel` grid spanning
`[-maximum(v), maximum(v)]` with an exact edge at `v_parallel = 0`.
"""
function default_vpar_edges(v::AbstractVector; n_edges::Int = 75)
    vmax = maximum(v)
    positive_edges = collect(range(0.0, vmax; length = n_edges))
    return vcat(-reverse(positive_edges[2:end]), positive_edges)
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

"""
    make_psd_from_AURORA(sim_dir; compute=:f_only, vpar_edges=nothing)

Load the AURORA simulation result from `sim_dir` and compute full phase-space
density `f`, the reduced distribution `F(v_parallel)`, or both.

If `vpar_edges` is `nothing` and `F` is requested, the function uses an
automatic symmetric uniform interval grid.
"""
function make_psd_from_AURORA(
    sim_dir::AbstractString;
    compute::Symbol = :f_only,
    vpar_edges::Union{Nothing, AbstractVector} = nothing,
)
    data = load_results(sim_dir)
    grids = psd_grids(data.E_centers, data.ΔE, data.μ_lims)

    if compute ∉ (:f_only, :F_only, :both)
        throw(ArgumentError("compute must be one of :f_only, :F_only, or :both"))
    end

    f_result =
        compute in (:f_only, :both) ?
        (f = compute_f(data.Ie, grids.ΔE_J, grids.ΔΩ, grids.v),) :
        NamedTuple()

    F_result =
        compute in (:F_only, :both) ?
        compute_F(data.Ie, data.μ_lims, grids.v; vpar_edges = vpar_edges) :
        NamedTuple()

    return merge(
        f_result,
        F_result,
        (
            v = grids.v,
            v_par = grids.v_par,
            v_perp = grids.v_perp,
            ΔE_J = grids.ΔE_J,
            BeamWeight = grids.ΔΩ,
            μ_center = grids.μ_center,
            E_edges = data.E_edges,
            E_centers = data.E_centers,
            t_run = data.t,
            h_atm = data.h_atm,
            μ_lims = data.μ_lims,
        ),
    )
end


# ======================================================================================== #
#                                  I/O AND WRAPPERS                                      #
# ======================================================================================== #

"""
    make_psd_file(directory_to_process; compute=:both, vpar_edges=nothing)

Read `simulation_data.nc` from `directory_to_process`, convert electron flux to
phase-space density, and write results to `analysis/psd.nc`.

# Keyword Arguments
- `compute`: one of `:f_only`, `:F_only`, or `:both`.
- `vpar_edges`: custom `v_parallel` bin edges [m/s]. If `nothing`, an automatic
  symmetric uniform grid is used spanning `[-maximum(v), maximum(v)]`.
"""
function make_psd_file(
    directory_to_process;
    compute::Symbol = :both,
    vpar_edges::Union{Nothing, AbstractVector} = nothing,
)
    println("Converting Ie to PSD.")

    result = make_psd_from_AURORA(directory_to_process; compute, vpar_edges)

    analysis_dir = joinpath(directory_to_process, "analysis")
    mkpath(analysis_dir)
    savefile = joinpath(analysis_dir, "psd.nc")
    write_psd_result(savefile, result)

    return nothing
end


function write_psd_result(savefile::AbstractString, res)
    h_atm = res.h_atm
    t_run = res.t_run
    Nz    = length(h_atm)
    Nt    = length(t_run)
    nE    = length(res.E_centers)
    nμ    = length(res.μ_lims) - 1

    NCDataset(savefile, "c") do ds
        defDim(ds, "altitude",          Nz)
        defDim(ds, "time",              Nt)
        defDim(ds, "energy",            nE)
        defDim(ds, "pitch_angle",       nμ)
        defDim(ds, "pitch_angle_bounds", nμ + 1)
        defDim(ds, "energy_bounds",      nE + 1)

        alt_v = defVar(ds, "altitude", Float64, ("altitude",);
                       attrib=["units" => "m", "long_name" => "altitude"])
        alt_v[:] = h_atm

        t_v = defVar(ds, "time", Float64, ("time",);
                     attrib=["units" => "s", "long_name" => "simulation time"])
        t_v[:] = t_run

        en_v = defVar(ds, "energy", Float64, ("energy",);
                      attrib=["units" => "eV", "long_name" => "energy bin center"])
        en_v[:] = res.E_centers

        ee_v = defVar(ds, "energy_edges", Float64, ("energy_bounds",);
                      attrib=["units" => "eV", "long_name" => "energy bin edges"])
        ee_v[:] = res.E_edges

        pa_v = defVar(ds, "pitch_angle", Float64, ("pitch_angle",);
                      attrib=["units" => "1", "long_name" => "cosine of pitch angle (beam center)"])
        pa_v[:] = res.μ_center

        ml_v = defVar(ds, "mu_lims", Float64, ("pitch_angle_bounds",);
                      attrib=["units" => "1", "long_name" => "pitch-angle cosine bin boundaries"])
        ml_v[:] = res.μ_lims

        v_v = defVar(ds, "v", Float64, ("energy",);
                     attrib=["units" => "m s-1", "long_name" => "electron speed"])
        v_v[:] = res.v

        bw_v = defVar(ds, "beam_weight", Float64, ("pitch_angle",);
                      attrib=["units" => "sr", "long_name" => "solid-angle beam weight"])
        bw_v[:] = res.BeamWeight

        dE_v = defVar(ds, "dE_J", Float64, ("energy",);
                      attrib=["units" => "J", "long_name" => "energy bin width in Joules"])
        dE_v[:] = res.ΔE_J

        if hasproperty(res, :f)
            f_v = defVar(ds, "f", Float32, ("altitude", "pitch_angle", "time", "energy");
                         deflatelevel=4,
                         chunksizes=(Nz, nμ, 1, nE),
                         attrib=["units"     => "s3 m-6",
                                  "long_name" => "phase-space density"])
            f_v[:, :, :, :] = Float32.(res.f)
        end

        if hasproperty(res, :F)
            Nvpar = length(res.vpar_centers)
            defDim(ds, "vpar",       Nvpar)
            defDim(ds, "vpar_bounds", Nvpar + 1)

            vpc_v = defVar(ds, "vpar_centers", Float64, ("vpar",);
                           attrib=["units" => "m s-1", "long_name" => "v_parallel bin centers"])
            vpc_v[:] = res.vpar_centers

            vpe_v = defVar(ds, "vpar_edges", Float64, ("vpar_bounds",);
                           attrib=["units" => "m s-1", "long_name" => "v_parallel bin edges"])
            vpe_v[:] = res.vpar_edges

            dv_v = defVar(ds, "dvpar", Float64, ("vpar",);
                          attrib=["units" => "m s-1", "long_name" => "v_parallel bin widths"])
            dv_v[:] = res.Δvpar

            F_v = defVar(ds, "F", Float32, ("vpar", "altitude", "time");
                         deflatelevel=4,
                         attrib=["units"     => "s m-4",
                                  "long_name" => "reduced distribution function F(v_parallel)"])
            F_v[:, :, :] = Float32.(res.F)
        end
    end
end

"""
    make_psd_file(sim::AuroraSimulation; kwargs...)

Convenience wrapper that calls [`make_psd_file`](@ref) on `sim.output.savedir`.
"""
make_psd_file(sim::AuroraSimulation; kwargs...) = make_psd_file(sim.output.savedir; kwargs...)
