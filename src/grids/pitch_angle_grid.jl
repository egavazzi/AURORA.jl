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
    θ_lims_vec = collect(θ_lims)

    maximum(θ_lims_vec) == 180 || throw(ArgumentError(
        "θ_lims must include 180° (field-aligned downward). " *
        "Got maximum of $(maximum(θ_lims_vec))°. " *
        "Example of valid input: 180:-10:0."
    ))
    minimum(θ_lims_vec) == 0 || throw(ArgumentError(
        "θ_lims must include 0° (field-aligned upward). " *
        "Got minimum of $(minimum(θ_lims_vec))°. " *
        "Example of valid input: 180:-10:0."
    ))
    issorted(θ_lims_vec, rev=true) || throw(ArgumentError(
        "θ_lims must be in descending order (e.g., 180:-10:0). " *
        "Got: $(θ_lims_vec)"
    ))

    μ_lims = cosd.(θ_lims_vec)
    μ_center = mu_avg(θ_lims_vec)
    FT = eltype(μ_lims)
    return PitchAngleGrid{FT, Vector{FT}}(
        collect(FT, θ_lims_vec), μ_lims, μ_center, length(θ_lims_vec) - 1
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
