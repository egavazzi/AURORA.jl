using SparseArrays: SparseArrays, spdiagm
using KLU: KLU, klu

mutable struct SolverCache
    KLU::KLU.KLUFactorization{Float64, Int64}
    initialized::Bool
    Mlhs::SparseArrays.SparseMatrixCSC{Float64, Int64}
    Mrhs::SparseArrays.SparseMatrixCSC{Float64, Int64}
    indices_lhs::Matrix{BlockIndices}     # nzval index map for Mlhs  (SS + CN)
    indices_rhs::Matrix{BlockIndices}     # nzval index map for Mrhs  (CN only)
    op_diags::OperatorDiagonals           # dense diagonals of Ddz_Up, Ddz_Down, Ddiffusion
    rhs::Vector{Float64}                  # reusable RHS buffer for Crank-Nicolson time-stepping
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
    rhs = Float64[]

    return SolverCache(KLU_factorization, false, Mlhs, Mrhs, indices_lhs, indices_rhs,
                       op_diags, rhs)
end

mutable struct DegradationCache{N}
    ionization_source_sum::Matrix{Float64}           # total incident flux summed over beams (shape: n_z x n_t)
    thermal_e_loss::Vector{Float64}
    Ie_scatter::Matrix{Float64}
    secondary_e_flux::NTuple{N, Matrix{Float64}}     # isotropic secondary e- flux per species (shape: n_z·n_μ x n_t)
    primary_e_flux::NTuple{N, Matrix{Float64}}       # forward primary e- flux per species    (shape: n_z·n_μ x n_t)
    secondary_e_spectrum::NTuple{N, Vector{Float64}} # energy spectrum weighting for secondaries per species (shape: n_E)
    primary_e_spectrum::NTuple{N, Vector{Float64}}   # energy spectrum weighting for primaries per species   (shape: n_E)
end


function DegradationCache{N}(n_μ::Int, n_t::Int, n_z::Int, n_E::Int) where {N}
    ionization_source_sum = Matrix{Float64}(undef, n_z, n_t)
    thermal_e_loss = Vector{Float64}(undef, n_z)
    Ie_scatter = Matrix{Float64}(undef, n_z * n_μ, n_t)

    secondary_e_flux = ntuple(_ -> zeros(n_z * n_μ, n_t), Val(N))
    primary_e_flux = ntuple(_ -> zeros(n_z * n_μ, n_t), Val(N))
    secondary_e_spectrum = ntuple(_ -> zeros(n_E), Val(N))
    primary_e_spectrum = ntuple(_ -> zeros(n_E), Val(N))

    return DegradationCache(ionization_source_sum, thermal_e_loss, Ie_scatter,
                            secondary_e_flux, primary_e_flux,
                            secondary_e_spectrum, primary_e_spectrum)
end

struct SimulationCache{D<:DegradationCache, TL, B}
    solver::SolverCache
    degradation::D
    matrices::TransportMatrices
    Ie::Array{Float64, 3}
    Ie_save::Array{Float64, 3}
    I0::Matrix{Float64}
    Ie_top::Array{Float64, 3}
    t_loop::TL
    B2B_fragment::B
end
