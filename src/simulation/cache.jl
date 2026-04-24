using SparseArrays: SparseArrays, spdiagm
using KLU: KLU, klu

mutable struct SolverCache
    KLU::KLU.KLUFactorization{Float64, Int64}
    Mlhs::SparseArrays.SparseMatrixCSC{Float64, Int64}
    Mrhs::SparseArrays.SparseMatrixCSC{Float64, Int64}
    indices_lhs::Matrix{BlockIndices}     # nzval index map for Mlhs  (SS + CN)
    indices_rhs::Matrix{BlockIndices}     # nzval index map for Mrhs  (CN only)
    op_diags::OperatorDiagonals           # dense diagonals of Ddz_Up, Ddz_Down, Ddiffusion
end

function SolverCache()
    KLU_factorization = klu(spdiagm(0 => ones(1)))
    Mlhs = spdiagm(0 => ones(1))
    Mrhs = spdiagm(0 => ones(1))
    dummy_bi = BlockIndices(Int[], Int[], Int[], 0, 0, 0)
    indices_lhs = fill(dummy_bi, 1, 1)
    indices_rhs = fill(dummy_bi, 1, 1)
    op_diags = OperatorDiagonals(
        [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0])

    return SolverCache(KLU_factorization, Mlhs, Mrhs, indices_lhs, indices_rhs, op_diags)
end

mutable struct DegradationCache
    n_repeated_over_μt::Vector{Matrix{Float64}}
    n_repeated_over_t::Vector{Matrix{Float64}}
    Ie_scatter::Matrix{Float64}
    secondary_e_flux::Vector{Matrix{Float64}}     # isotropic secondary e- flux per species (shape: n_z·n_μ x n_t)
    primary_e_flux::Vector{Matrix{Float64}}       # forward primary e- flux per species    (shape: n_z·n_μ x n_t)
    secondary_e_spectrum::Vector{Vector{Float64}} # energy spectrum weighting for secondaries per species (shape: n_E)
    primary_e_spectrum::Vector{Vector{Float64}}   # energy spectrum weighting for primaries per species   (shape: n_E)
end

function DegradationCache()
    n_repeated_over_μt = [zeros(1, 1) for _ in 1:3]
    n_repeated_over_t = [zeros(1, 1) for _ in 1:3]
    Ie_scatter = zeros(1, 1)
    secondary_e_flux     = [zeros(1, 1) for _ in 1:3]
    primary_e_flux       = [zeros(1, 1) for _ in 1:3]
    secondary_e_spectrum = [zeros(1) for _ in 1:3]
    primary_e_spectrum   = [zeros(1) for _ in 1:3]
    return DegradationCache(n_repeated_over_μt, n_repeated_over_t, Ie_scatter,
                            secondary_e_flux, primary_e_flux,
                            secondary_e_spectrum, primary_e_spectrum)
end

function DegradationCache(n_neutrals, n_μ::Int, n_t::Int, n_z::Int, n_E::Int)
    number_species = length(n_neutrals)

    n_repeated_over_μt = Vector{Matrix{Float64}}(undef, number_species)
    n_repeated_over_t = Vector{Matrix{Float64}}(undef, number_species)

    for i_species in 1:number_species
        n = n_neutrals[i_species]

        n_repeated_over_μt[i_species] = Matrix{Float64}(undef, n_z * n_μ, n_t)
        for i_t in 1:n_t
            for i_μ in 1:n_μ
                @views n_repeated_over_μt[i_species][(i_μ - 1) * n_z .+ (1:n_z), i_t] .= n
            end
        end

        n_repeated_over_t[i_species] = Matrix{Float64}(undef, n_z, n_t)
        for i_t in 1:n_t
            @views n_repeated_over_t[i_species][:, i_t] .= n
        end
    end

    Ie_scatter = Matrix{Float64}(undef, n_z * n_μ, n_t)

    secondary_e_flux     = Vector{Matrix{Float64}}(undef, number_species)
    primary_e_flux       = Vector{Matrix{Float64}}(undef, number_species)
    secondary_e_spectrum = Vector{Vector{Float64}}(undef, number_species)
    primary_e_spectrum   = Vector{Vector{Float64}}(undef, number_species)

    for i_species in 1:number_species
        secondary_e_flux[i_species]     = zeros(n_z * n_μ, n_t)
        primary_e_flux[i_species]       = zeros(n_z * n_μ, n_t)
        secondary_e_spectrum[i_species] = zeros(n_E)
        primary_e_spectrum[i_species]   = zeros(n_E)
    end

    return DegradationCache(n_repeated_over_μt, n_repeated_over_t, Ie_scatter,
                            secondary_e_flux, primary_e_flux,
                            secondary_e_spectrum, primary_e_spectrum)
end

struct SimulationCache
    solver::SolverCache
    degradation::DegradationCache
    cascading::CascadingCache
    matrices::TransportMatrices
    Ie::Array{Float64, 3}
    Ie_save::Array{Float64, 3}
    I0::Matrix{Float64}
    Ie_top::Array{Float64, 3}
    t_loop
    phase_fcn_neutrals
    B2B_fragment
end
