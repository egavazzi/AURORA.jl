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
end

"""
    TransportMatrices(n_altitude, n_angle, n_time, n_energy)

Construct an empty TransportMatrices container with zeros.
"""
function TransportMatrices(n_altitude::Int, n_angle::Int, n_time::Int, n_energy::Int)
    A = zeros(Float64, n_altitude)
    B = zeros(Float64, n_altitude, n_angle, n_angle)
    D = zeros(Float64, n_energy, n_angle)
    Q = zeros(Float64, n_altitude * n_angle, n_time, n_energy)
    Ddiffusion = SparseArrays.spzeros(Float64, n_altitude, n_altitude)
    return TransportMatrices(A, B, D, Q, Ddiffusion)
end
