"""
    PitchAngleGrid{FT, V<:AbstractVector{FT}} <: AbstractGrid

Pitch-angle grid for electron beams.
"""
struct PitchAngleGrid{FT, V<:AbstractVector{FT}} <: AbstractGrid
    θ_lims::V
    μ_lims::V
    μ_center::V
    n_beams::Int
end

function PitchAngleGrid(θ_lims)
    validate_θ_lims(θ_lims)
    μ_lims = cosd.(collect(θ_lims))
    μ_center = mu_avg(θ_lims)
    FT = eltype(μ_lims)
    return PitchAngleGrid{FT, Vector{FT}}(
        collect(FT, θ_lims), μ_lims, μ_center, length(θ_lims) - 1
    )
end

function Base.show(io::IO, grid::PitchAngleGrid)
    print(io, "PitchAngleGrid($(grid.n_beams) beams)")
end

function Base.show(io::IO, ::MIME"text/plain", grid::PitchAngleGrid)
    println(io, "PitchAngleGrid:")
    println(io, "├── Beams: $(grid.n_beams)")
    print(io, "└── θ_lims: $(grid.θ_lims[1])° to $(grid.θ_lims[end])°")
end
