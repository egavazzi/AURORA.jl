@testitem "calculate_volume_excitation basic Q = Ie x σ x n" begin
    using AURORA

    n_z = 3
    n_t = 2
    n_E = 4

    z = collect(1.0:n_z) .* 1e5  # altitudes in m
    t = collect(1.0:n_t)

    # Uniform flux of 1.0 everywhere
    Ie_ztE_omni = ones(n_z, n_t, n_E)

    # Cross-section: linearly increasing with energy index
    σ = collect(1.0:n_E)

    # Neutral density: different at each altitude
    n = collect(1.0:n_z) .* 1e10

    Q = AURORA.calculate_volume_excitation(z, t, Ie_ztE_omni, σ, n)

    @test size(Q) == (n_z, n_t)
    # Q[iz, it] = sum_E(Ie[iz, it, :] * σ[:]) * n[iz]
    # With Ie = 1 everywhere: Q[iz, it] = sum(σ) * n[iz]
    expected_sum_σ = sum(σ)
    for iz in 1:n_z, it in 1:n_t
        @test Q[iz, it] ≈ expected_sum_σ * n[iz]
    end
end

@testitem "q2colem steady-state integrates correctly" begin
    using AURORA

    # Constant Q = 1.0 at all altitudes → integral should equal the total height span
    z = collect(range(100e3, 200e3; length = 50))
    Q = ones(length(z))
    t_scalar = 0.0

    I = AURORA.q2colem(t_scalar, z, Q)

    @test length(I) == 1
    # Trapezoidal integration of constant 1.0 over [100e3, 200e3] = 100e3
    @test I[1] ≈ 100e3 rtol = 1e-10
end

@testitem "q2colem time-dependent returns correct shape" begin
    using AURORA

    n_z = 20
    n_t = 10
    z = collect(range(100e3, 200e3; length = n_z))
    t = collect(range(0.0, 1.0; length = n_t))

    # Constant Q everywhere → no time-delay effect matters much
    Q = ones(n_z, n_t)

    I = AURORA.q2colem(t, z, Q)

    @test length(I) == n_t
    @test all(I .>= 0)
    # For constant Q, the last time steps (where photon delays don't truncate the integral)
    # should be the same as steady state (integral of 1.0 over z = 100e3)
    @test I[end] ≈ 100e3 rtol = 0.01
end
