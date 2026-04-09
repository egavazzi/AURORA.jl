@testitem "calculate_heating_rate output shape and positivity" begin
    using AURORA

    n_z = 5
    n_t = 3
    n_E = 4

    z = collect(range(100e3, 300e3; length = n_z))
    t = collect(range(0.0, 1.0; length = n_t))
    E_centers = collect(range(10.0, 500.0; length = n_E))

    # Uniform positive flux
    Ie_ztE_omni = ones(n_z, n_t, n_E) .* 1e8

    # Typical ionospheric values
    ne = fill(1e11, n_z)   # m⁻³
    Te = fill(2000.0, n_z) # K

    hr = AURORA.calculate_heating_rate(z, t, Ie_ztE_omni, E_centers, ne, Te)

    @test size(hr) == (n_z, n_t)
    # Heating rate should be positive when there is positive flux flowing
    # through a plasma with thermal electrons
    @test all(hr .> 0)
end
