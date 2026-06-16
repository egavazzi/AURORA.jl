using SparseArrays: SparseArrays

"""
    TransportMatrices

Container for the matrices used in the electron transport equations.
"""
mutable struct TransportMatrices
    A::Vector{Float64}
    B::Array{Float64, 3}
    D::Array{Float64, 2}
    Q::Array{Float64, 3}
    Ddiffusion::SparseArrays.SparseMatrixCSC{Float64, Int64}

    Le::Vector{Float64}                              # thermal e- energy loss (length n_z)
    B2B_elastic::Matrix{Float64}                     # elastic beam-to-beam (n_angle × n_angle)
    B2B_inelastic_neutrals::Vector{Matrix{Float64}}  # inelastic beam-to-beam per species
    phase_fcn_e::Vector{Float64}                     # 3D elastic phase function scratch
    phase_fcn_i::Vector{Float64}                     # 3D inelastic phase function scratch
end

"""
    TransportMatrices(n_altitude, n_angle, n_time, n_energy)

Construct an empty TransportMatrices container with zeros.

The scratch buffers that depend on the number of species and on the scattering grid
(`B2B_inelastic_neutrals`, `phase_fcn_e`, `phase_fcn_i`) are created empty here and
sized later in [`initialize_transport_matrices`](@ref).
"""
function TransportMatrices(n_altitude::Int, n_angle::Int, n_time::Int, n_energy::Int)
    A = zeros(Float64, n_altitude)
    B = zeros(Float64, n_altitude, n_angle, n_angle)
    D = zeros(Float64, n_energy, n_angle)
    Q = zeros(Float64, n_altitude * n_angle, n_time, n_energy)
    Ddiffusion = SparseArrays.spzeros(Float64, n_altitude, n_altitude)

    Le = zeros(Float64, n_altitude)
    B2B_elastic = zeros(Float64, n_angle, n_angle)
    B2B_inelastic_neutrals = Matrix{Float64}[]
    phase_fcn_e = Float64[]
    phase_fcn_i = Float64[]

    return TransportMatrices(A, B, D, Q, Ddiffusion,
                             Le, B2B_elastic, B2B_inelastic_neutrals,
                             phase_fcn_e, phase_fcn_i)
end
