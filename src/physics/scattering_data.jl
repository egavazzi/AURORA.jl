"""
    ScatteringData{FT, A<:AbstractArray{FT}, M<:AbstractMatrix{FT}, V<:AbstractVector{FT}}

Pre-computed scattering probability matrices for beam-to-beam scattering.
"""
struct ScatteringData{FT, A<:AbstractArray{FT}, M<:AbstractMatrix{FT}, V<:AbstractVector{FT}}
    Pmu2mup::A
    BeamWeight_relative::M
    BeamWeight::V
    theta1::V
end

function ScatteringData(θ_lims; n_direction=720)
    validate_θ_lims(θ_lims)
    BeamWeight = beam_weight(θ_lims)
    Pmu2mup, _, BeamWeight_relative, θ1 = find_scattering_matrices(θ_lims, n_direction)
    FT = eltype(Pmu2mup)
    return ScatteringData{FT, typeof(Pmu2mup), typeof(BeamWeight_relative), typeof(BeamWeight)}(
        Pmu2mup, BeamWeight_relative, BeamWeight, vec(θ1)
    )
end

function Base.show(io::IO, sd::ScatteringData)
    n_beams = length(sd.BeamWeight)
    n_dir = length(sd.theta1)
    print(io, "ScatteringData($(n_beams) beams, $(n_dir) directions)")
end

function Base.show(io::IO, ::MIME"text/plain", sd::ScatteringData)
    n_beams = length(sd.BeamWeight)
    n_dir = length(sd.theta1)
    println(io, "ScatteringData:")
    println(io, "├── Beams: $(n_beams)")
    println(io, "├── Directions: $(n_dir)")
    print(io, "└── Pmu2mup: $(join(size(sd.Pmu2mup), 'x'))")
end
