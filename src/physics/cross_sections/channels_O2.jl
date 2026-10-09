"""
    default_channels_O2() → Vector{CollisionChannel}

Built-in O₂ collision channels.
"""
function default_channels_O2()
    return CollisionChannel[
        CollisionChannel("_OO3S", e_O2_OO3S, 15.6, 0),
        CollisionChannel("_9p97", e_O2_9p97, 9.97, 0),
        CollisionChannel("_8p4", e_O2_8p4, 8.4, 0),
        CollisionChannel("_6", e_O2_6, 6.0, 0),
        CollisionChannel("_4p5", e_O2_4p5, 4.5, 0),
        CollisionChannel("b1Sgp", e_O2b1Sgp, 1.627, 0),
        CollisionChannel("a1Dg", e_O2a1Dg, 0.977, 0),
        CollisionChannel("vib", e_O2vib, 0.193, 0),
        CollisionChannel("ionx2pg", e_O2ionx2pg, 12.072, 1),
        CollisionChannel("iona4pu", e_O2iona4pu, 16.1, 1),
        CollisionChannel("ion16p9", e_O2ion16p9, 16.9, 1),
        CollisionChannel("ionb4sgm", e_O2ionb4sgm, 18.2, 1),
        CollisionChannel("dion", e_O2dion, 18.9, 1),
        CollisionChannel("ddion", e_O2ddion, 32.51, 2),
    ]
end

"""
    default_elastic_cross_section_O2()

The built-in O₂ elastic cross section, `e_O2elastic`.
"""
default_elastic_cross_section_O2() = e_O2elastic

"""
    default_secondary_law_O2()

The built-in O₂ secondary-electron energy distribution `(E_s, E_p) -> 1 / (15.2² + E_s²)`,
for building the cascading transfer matrices.
"""
default_secondary_law_O2() = @law (E_s, E_p) -> 1.0 / (15.2^2 + E_s^2)
