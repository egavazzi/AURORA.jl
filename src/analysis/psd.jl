using MAT: matopen, matread
using LoopVectorization

const m_e = 9.10938356e-31
const e_charge = 1.602176634e-19

"""
    load_Ie(path)

Load a `IeFlickering-NN.mat` result file produced by AURORA.

Returns a named tuple with fields:
- `Ie`      : `[Nz, n_mu, Nt, nE]` number flux [m^-2 s^-1], integrated over ΔE and ΔΩ
- `E_edges` : energy bin edges [eV], length `nE + 1`
- `E_centers`: energy bin centers [eV], length `nE`
- `dE`      : energy bin widths [eV], length `nE`
- `μ_lims`  : cosine pitch-angle bin limits, length `n_mu + 1`
- `t_run`   : time vector [s], length `Nt`
- `h_atm`   : altitude grid [m], length `Nz`
"""
function load_Ie(path::AbstractString)
    file = matopen(path)
    Ie_raw = read(file, "Ie_ztE")
    E_edges = vec(read(file, "E_edges"))
    E_centers = vec(read(file, "E_centers"))
    ΔE = vec(read(file, "dE"))
    μ_lims = vec(read(file, "mu_lims"))
    t_run = vec(read(file, "t_run"))
    h_atm = vec(read(file, "h_atm"))
    close(file)

    nE = size(Ie_raw, 3)
    Nz = length(h_atm)
    nμ = length(μ_lims) - 1

    Ie = reshape(Ie_raw, Nz, nμ, length(t_run), nE)

    return (Ie = Ie, E_edges = E_edges, E_centers = E_centers, ΔE = ΔE, μ_lims = μ_lims, t_run = t_run, h_atm = h_atm)
end

"""
    psd_grids(E_centers, ΔE, μ_lims_cosine)

Compute velocity and pitch-angle grids needed by `compute_f` and `compute_F`.

- `E_centers` is the energy bin centers [eV] (length `nE`).
- `ΔE` is the energy bin widths [eV] (length `nE`).
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
    # Use a [Nz, Nvpar, Nt] intermediate layout for the computation.
    # This ensures iz is the innermost dimension for both arrays, which enables much more
    # efficient SIMD optimizations from @tturbo compared to the original [Nvpar, Nz, Nt] layout.
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

    # Permute [Nz, Nvpar, Nt] → [Nvpar, Nz, Nt]
    F = permutedims(F_local, (2, 1, 3))

    return (F = F, vpar_edges = vpar_edges, vpar_centers = vpar_centers, Δvpar = Δvpar)
end

"""
    make_psd_from_AURORA(path_to_file; compute=:f_only, vpar_edges=nothing)

Load an AURORA result file and compute full phase-space density `f`, the reduced
distribution `F(v_parallel)`, or both.

If `vpar_edges` is `nothing` and `F` is requested, the function uses an
automatic symmetric uniform interval grid.
"""
function make_psd_from_AURORA(
    path_to_file::AbstractString;
    compute::Symbol = :f_only,
    vpar_edges::Union{Nothing, AbstractVector} = nothing,
)
    data = load_Ie(path_to_file)
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
            t_run = data.t_run,
            h_atm = data.h_atm,
            μ_lims = data.μ_lims,
        ),
    )
end


# ======================================================================================== #
#                                  I/O AND WRAPPERS                                      #
# ======================================================================================== #

"""
    make_psd_file(directory_to_process; compute=:both, vpar_edges=nothing, output_prefix="psd")

Reads into a folder `directory_to_process` containing results from an AURORA.jl simulation,
loads particle flux files `IeFlickering-NN.mat`, converts them into phase-space density, and
writes one output file per input as `psd/psd-NN.mat`.

# Calling
`make_psd_file(directory_to_process; compute=:both, vpar_edges=nothing, output_prefix="psd")`

# Inputs
- `directory_to_process`: absolute or relative path to the simulation directory to process.

# Keyword Arguments
- `compute`: one of `:f_only`, `:F_only`, or `:both`.
- `vpar_edges`: custom `v_parallel` bin edges [m/s] used when computing `F`.
    If `nothing`, an automatic symmetric uniform interval grid is used, spanning
    `[-maximum(v), maximum(v)]` with an edge at `v_parallel = 0`.
- `output_prefix`: output filename prefix used as `<output_prefix>-NN.mat`.
"""
function make_psd_file(
    directory_to_process;
    compute::Symbol = :both,
    vpar_edges::Union{Nothing, AbstractVector} = nothing,
    output_prefix::AbstractString = "psd",
)
    println("Converting Ie to PSD.")
    files = readdir(directory_to_process, join = true)
    files_to_process = files[contains.(files, r"IeFlickering\-[0-9]+\.mat")]
    sort!(
        files_to_process,
        by = x -> parse(Int, match(r"IeFlickering-(\d+)\.mat", basename(x))[1]),
    )

    if isempty(files_to_process)
        @warn "No simulation results found in $directory_to_process. Skipping phase-space density calculations."
        return nothing
    end

    n_files = length(files_to_process)
    psd_dir = joinpath(directory_to_process, "psd")
    mkpath(psd_dir)

    for (i_file, file) in enumerate(files_to_process)
        match_result = match(r"IeFlickering-(\d+)\.mat", basename(file))
        @assert !isnothing(match_result) "Unexpected input filename: $file"
        file_id = match_result[1]

        print("\rProcessing PSD: $(i_file)/$(n_files)")
        flush(stdout)

        result = make_psd_from_AURORA(file; compute = compute, vpar_edges = vpar_edges)
        outfile = joinpath(psd_dir, "$(output_prefix)-$(file_id).mat")
        write_psd_result(outfile, result)
    end

    println()

    return nothing
end


function write_psd_result(outfile::AbstractString, res)
    matopen(outfile, "w") do io
        if hasproperty(res, :f)
            write(io, "f", res.f)
        end

        if hasproperty(res, :F)
            write(io, "F", res.F)
            write(io, "vpar_edges", res.vpar_edges)
            write(io, "vpar_centers", res.vpar_centers)
            write(io, "dvpar", res.Δvpar)
        end

        write(io, "v", res.v)
        write(io, "v_par", res.v_par)
        write(io, "v_perp", res.v_perp)
        write(io, "dE_J", res.ΔE_J)
        write(io, "BeamWeight", res.BeamWeight)
        write(io, "mu_center", res.μ_center)
        write(io, "E_edges", res.E_edges)
        write(io, "E_centers", res.E_centers)
        write(io, "t_run", res.t_run)
        write(io, "h_atm", res.h_atm)
        write(io, "mu_lims", res.μ_lims)
    end
end

"""
    make_psd_file(sim::AuroraSimulation; kwargs...)

Convenience wrapper that calls [`make_psd_file`](@ref) on `sim.savedir`.
"""
make_psd_file(sim::AuroraSimulation; kwargs...) = make_psd_file(sim.savedir; kwargs...)
