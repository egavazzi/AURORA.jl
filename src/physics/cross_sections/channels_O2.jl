"""
    default_collisions_O2() → NamedTuple

Built-in O₂ collision data `(; elastic_cross_section, channels, secondary_law)`. The
secondary-electron energy distribution is `(E_s, E_p) -> 1 / (15.2² + E_s²)`.
"""
function default_collisions_O2()
    elastic_cross_section = e_O2elastic
    channels = [
        CollisionChannel("_OO3S", e_O2_OO3S, 15.6, 0; source = "built-in"),
        CollisionChannel("_9p97", e_O2_9p97, 9.97, 0; source = "built-in"),
        CollisionChannel("_8p4", e_O2_8p4, 8.4, 0; source = "built-in"),
        CollisionChannel("_6", e_O2_6, 6.0, 0; source = "built-in"),
        CollisionChannel("_4p5", e_O2_4p5, 4.5, 0; source = "built-in"),
        CollisionChannel("b1Sgp", e_O2b1Sgp, 1.627, 0; source = "built-in"),
        CollisionChannel("a1Dg", e_O2a1Dg, 0.977, 0; source = "built-in"),
        CollisionChannel("vib", e_O2vib, 0.193, 0; source = "built-in"),
        CollisionChannel("ionx2pg", e_O2ionx2pg, 12.072, 1; source = "built-in"),
        CollisionChannel("iona4pu", e_O2iona4pu, 16.1, 1; source = "built-in"),
        CollisionChannel("ion16p9", e_O2ion16p9, 16.9, 1; source = "built-in"),
        CollisionChannel("ionb4sgm", e_O2ionb4sgm, 18.2, 1; source = "built-in"),
        CollisionChannel("dion", e_O2dion, 18.9, 1; source = "built-in"),
        CollisionChannel("ddion", e_O2ddion, 32.51, 2; source = "built-in"),
    ]
    secondary_law = @law (E_s, E_p) -> 1.0 / (15.2^2 + E_s^2)
    return (; elastic_cross_section, channels, secondary_law)
end
