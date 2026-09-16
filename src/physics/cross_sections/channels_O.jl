"""
    default_channels_O() → Vector{CollisionChannel}

The built-in atomic-oxygen inelastic collision channels, in the row order of the
`cross_sections` and `excitation_levels` matrices of an O species.

The ground-state fine-structure transitions (³P₂→³P₀ at 0.0281 eV, ³P₂→³P₁ at 0.0196 eV and
³P₁→³P₀ at 0.0085 eV) are absent: no cross section is available for them.
"""
function default_channels_O()
    return CollisionChannel[
        CollisionChannel("1D", e_O1D, 1.967, 0),
        CollisionChannel("1S", e_O1S, 4.19, 0),
        CollisionChannel("3s5S0", e_O3s5S0, 9.14, 0),
        CollisionChannel("3s3S0", e_O3s3S0, 9.521, 0),
        CollisionChannel("3p5P", e_O3p5P, 10.74, 0),
        CollisionChannel("3sp3D0", e_O3sp3D0, 12.54, 0),
        CollisionChannel("3p3P", e_O3p3P, 10.99, 0),
        CollisionChannel("ion4S0", e_Oion4S0, 13.618, 1),
        CollisionChannel("ion2D0", e_Oion2D0, 16.9, 1),
        CollisionChannel("ion2P0", e_Oion2P0, 18.6, 1),
        CollisionChannel("ionion", e_Oionion, 28.5, 2),
    ]
end

"""
    default_elastic_cross_section_O()

The built-in atomic-oxygen elastic cross section, [`e_Oelastic`](@ref).
"""
default_elastic_cross_section_O() = e_Oelastic

"""
    default_secondary_law_O()

The built-in atomic-oxygen secondary-electron energy distribution, an [`OSecondaryLaw`](@ref)
interpolating its parameters in the primary energy.
"""
function default_secondary_law_O()
    energy_params = [100.0, 200, 500, 1000, 2000]  # eV
    A_params = [12.6, 13.7, 14.1, 14.0, 13.7]
    B_params = [7.18, 4.97, 2.75, 1.69, 1.02] .* 1e-22
    return OSecondaryLaw(energy_params, A_params, B_params)
end
