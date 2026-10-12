"""
    default_collisions_N2() → NamedTuple

Built-in N₂ collision data `(; elastic_cross_section, channels, secondary_law)`. The
secondary-electron energy distribution is `(E_s, E_p) -> 1 / (11.4² + E_s²)`.
"""
function default_collisions_N2()
    elastic_cross_section = e_N2elastic
    channels = [
        CollisionChannel("rot0_2", e_N2rot0_2, 0.00148010556, 0; source = "built-in"),
        CollisionChannel("rot0_4", e_N2rot0_4, 0.004933884, 0; source = "built-in"),
        CollisionChannel("rot0_6", e_N2rot0_6, 0.01036181244, 0; source = "built-in"),
        CollisionChannel("rot0_8", e_N2rot0_8, 0.01776464064, 0; source = "built-in"),
        CollisionChannel("vib0_1", e_N2vib0_1, 0.2888, 0; source = "built-in"),
        CollisionChannel("vib0_2", e_N2vib0_2, 0.5742, 0; source = "built-in"),
        CollisionChannel("vib0_3", e_N2vib0_3, 0.8559, 0; source = "built-in"),
        CollisionChannel("vib0_4", e_N2vib0_4, 1.1342, 0; source = "built-in"),
        CollisionChannel("vib0_5", e_N2vib0_5, 1.4088, 0; source = "built-in"),
        CollisionChannel("vib0_6", e_N2vib0_6, 1.6801, 0; source = "built-in"),
        CollisionChannel("vib0_7", e_N2vib0_7, 1.9475, 0; source = "built-in"),
        CollisionChannel("a3sup", e_N2a3sup, 6.1688, 0; source = "built-in"),
        CollisionChannel("b3pg", e_N2b3pg, 7.3532, 0; source = "built-in"),
        CollisionChannel("w3du", e_N2w3du, 7.3622, 0; source = "built-in"),
        CollisionChannel("bp3sum", e_N2bp3sum, 8.1647, 0; source = "built-in"),
        CollisionChannel("ap1sum", e_N2ap1sum, 8.3987, 0; source = "built-in"),
        CollisionChannel("w1du", e_N2w1du, 8.8895, 0; source = "built-in"),
        CollisionChannel("e3sgp", e_N2e3sgp, 11.875, 0; source = "built-in"),
        CollisionChannel("ab1sgp", e_N2ab1sgp, 12.255, 0; source = "built-in"),
        CollisionChannel("a1pg", e_N2a1pg, 8.5489, 0; source = "built-in"),
        CollisionChannel("c3pu", e_N2c3pu, 11.032, 0; source = "built-in"),
        CollisionChannel("bp1sup", e_N2bp1sup, 12.85, 0; source = "built-in"),
        CollisionChannel("cp1sup", e_N2cp1sup, 12.94, 0; source = "built-in"),
        CollisionChannel("cp3pu", e_N2cp3pu, 12.08, 0; source = "built-in"),
        CollisionChannel("d3sup", e_N2d3sup, 12.85, 0; source = "built-in"),
        CollisionChannel("f3pu", e_N2f3pu, 12.75, 0; source = "built-in"),
        CollisionChannel("g3pu", e_N2g3pu, 12.8, 0; source = "built-in"),
        CollisionChannel("M1M2", e_N2M1M2, 13.15, 0; source = "built-in"),
        CollisionChannel("o1pu", e_N2o1pu, 13.1, 0; source = "built-in"),
        CollisionChannel("dissociation", e_N2dissociation, 20.6, 0; source = "built-in"),
        CollisionChannel("ionx2sgp", e_N2ionx2sgp, 15.581, 1; source = "built-in"),
        CollisionChannel("iona2pu", e_N2iona2pu, 16.73, 1; source = "built-in"),
        CollisionChannel("ionb2sup", e_N2ionb2sup, 18.75, 1; source = "built-in"),
        CollisionChannel("dion", e_N2dion, 24.0, 1; source = "built-in"),
        CollisionChannel("ddion", e_N2ddion, 42.0, 2; source = "built-in"),
    ]
    secondary_law = @law (E_s, E_p) -> 1.0 / (11.4^2 + E_s^2)
    return (; elastic_cross_section, channels, secondary_law)
end
