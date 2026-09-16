"""
    default_channels_N2() → Vector{CollisionChannel}

The built-in N₂ inelastic collision channels, in the row order of the `cross_sections` and
`excitation_levels` matrices of an N₂ species.
"""
function default_channels_N2()
    return CollisionChannel[
        CollisionChannel("rot0_2", e_N2rot0_2, 0.00148010556, 0),
        CollisionChannel("rot0_4", e_N2rot0_4, 0.004933884, 0),
        CollisionChannel("rot0_6", e_N2rot0_6, 0.01036181244, 0),
        CollisionChannel("rot0_8", e_N2rot0_8, 0.01776464064, 0),
        CollisionChannel("vib0_1", e_N2vib0_1, 0.2888, 0),
        CollisionChannel("vib0_2", e_N2vib0_2, 0.5742, 0),
        CollisionChannel("vib0_3", e_N2vib0_3, 0.8559, 0),
        CollisionChannel("vib0_4", e_N2vib0_4, 1.1342, 0),
        CollisionChannel("vib0_5", e_N2vib0_5, 1.4088, 0),
        CollisionChannel("vib0_6", e_N2vib0_6, 1.6801, 0),
        CollisionChannel("vib0_7", e_N2vib0_7, 1.9475, 0),
        CollisionChannel("a3sup", e_N2a3sup, 6.1688, 0),
        CollisionChannel("b3pg", e_N2b3pg, 7.3532, 0),
        CollisionChannel("w3du", e_N2w3du, 7.3622, 0),
        CollisionChannel("bp3sum", e_N2bp3sum, 8.1647, 0),
        CollisionChannel("ap1sum", e_N2ap1sum, 8.3987, 0),
        CollisionChannel("w1du", e_N2w1du, 8.8895, 0),
        CollisionChannel("e3sgp", e_N2e3sgp, 11.875, 0),
        CollisionChannel("ab1sgp", e_N2ab1sgp, 12.255, 0),
        CollisionChannel("a1pg", e_N2a1pg, 8.5489, 0),
        CollisionChannel("c3pu", e_N2c3pu, 11.032, 0),
        CollisionChannel("bp1sup", e_N2bp1sup, 12.85, 0),
        CollisionChannel("cp1sup", e_N2cp1sup, 12.94, 0),
        CollisionChannel("cp3pu", e_N2cp3pu, 12.08, 0),
        CollisionChannel("d3sup", e_N2d3sup, 12.85, 0),
        CollisionChannel("f3pu", e_N2f3pu, 12.75, 0),
        CollisionChannel("g3pu", e_N2g3pu, 12.8, 0),
        CollisionChannel("M1M2", e_N2M1M2, 13.15, 0),
        CollisionChannel("o1pu", e_N2o1pu, 13.1, 0),
        CollisionChannel("dissociation", e_N2dissociation, 20.6, 0),
        CollisionChannel("ionx2sgp", e_N2ionx2sgp, 15.581, 1),
        CollisionChannel("iona2pu", e_N2iona2pu, 16.73, 1),
        CollisionChannel("ionb2sup", e_N2ionb2sup, 18.75, 1),
        CollisionChannel("dion", e_N2dion, 24.0, 1),
        CollisionChannel("ddion", e_N2ddion, 42.0, 2),
    ]
end

"""
    default_elastic_cross_section_N2()

The built-in N₂ elastic cross section, [`e_N2elastic`](@ref).
"""
default_elastic_cross_section_N2() = e_N2elastic

"""
    default_secondary_law_N2()

The built-in N₂ secondary-electron energy distribution `(E_s, E_p) -> 1 / (11.4² + E_s²)`,
for building the cascading transfer matrices.
"""
default_secondary_law_N2() = @law (E_s, E_p) -> 1.0 / (11.4^2 + E_s^2)
