# Non-ionizing degradation: `add_inelastic_collisions!` places the electrons that leave an
# energy bin into the lower bins they arrive in, and `update_B!` retains the rest in place.
@testmodule DegradationSetup begin
    using AURORA

    # Uniform 1 eV grid with its lowest edge at 1 eV, so the sub-floor cases are easy to read.
    E_EDGES = collect(1.0:1.0:11.0)
    ΔE = diff(E_EDGES)
    E_CENTERS = E_EDGES[1:end-1] .+ ΔE ./ 2
    GRID = AURORA.EnergyGrid{Float64, Vector{Float64}}(E_EDGES, E_CENTERS, ΔE,
                                                       length(ΔE), E_EDGES[end])

    # One altitude, one beam, one time step, a single non-ionizing channel with unit cross
    # section, and a unit flux in bin `iE`. Returns the energy spectrum deposited in Q.
    function degrade(iE, E_loss)
        n_E = GRID.n
        Q = zeros(1, 1, n_E)
        Ie = zeros(1, 1, n_E)
        Ie[1, 1, iE] = 1.0
        E_levels = [0.0 0.0; E_loss 0.0]
        σ = zeros(2, n_E)
        σ[2, :] .= 1.0
        workspace = (; Ie_scatter = zeros(1, 1))
        AURORA.add_inelastic_collisions!(Q, Ie, [0.0], [1.0], σ, E_levels, ones(1, 1),
                                         ones(1, 1, 1), GRID, iE, workspace)
        return Q[1, 1, :]
    end

    # Fraction of the bin's electrons that leave it, which `update_B!` complements with
    # max(0, 1 - E_loss/ΔE).
    leaving_fraction(iE, E_loss) = min(1, E_loss / ΔE[iE])
end

@testitem "Inelastic degradation on grid removes the channel's energy loss" setup=[DegradationSetup] begin
    E_CENTERS = DegradationSetup.E_CENTERS

    for (iE, E_loss) in ((8, 2.5), (8, 0.4), (6, 1.0), (10, 3.0))
        placed = DegradationSetup.degrade(iE, E_loss)
        leaving = DegradationSetup.leaving_fraction(iE, E_loss)

        # Every leaving electron is placed: nothing is lost above the grid floor.
        @test sum(placed) ≈ leaving
        # And the energy they carry away is the channel's energy loss.
        mean_arrival = sum(placed .* E_CENTERS) / sum(placed)
        @test leaving * (E_CENTERS[iE] - mean_arrival) ≈ E_loss
        # Nothing is placed at or above the bin itself.
        @test all(placed[iE:end] .== 0)
    end
end

@testitem "Inelastic degradation below the grid floor places only the on-grid share" setup=[DegradationSetup] begin
    # Bin 2 spans [2, 3] eV; losing 1.5 eV maps the leaving electrons onto [0.5, 1.5] eV, of
    # which only [1, 1.5] is on the grid.
    placed = DegradationSetup.degrade(2, 1.5)
    leaving = DegradationSetup.leaving_fraction(2, 1.5)

    @test placed[1] ≈ leaving * 0.5
    @test sum(placed) ≈ leaving * 0.5
    @test sum(placed) < leaving
    @test all(placed[2:end] .== 0)

    # Electrons leaving the lowest bin have nowhere on the grid to go.
    @test all(DegradationSetup.degrade(1, 0.4) .== 0)
    @test all(DegradationSetup.degrade(1, 5.0) .== 0)
end

@testitem "Scaled beam-to-beam transfers conserve the cell-weighted column" begin
    # Non-uniform field-line grid and two beams in each direction.
    s = cumsum([0.0, 1.0, 1.5, 2.5, 4.0, 6.5, 10.0])
    μ = [-0.8, -0.3, 0.3, 0.8]
    n_z, n_μ = length(s), length(μ)
    F = AURORA.upwind_cell_ratio(s, μ)
    w = AURORA.column_weights(s)
    wrow(iz, i) = μ[i] < 0 ? w.down[iz] : w.up[iz]

    @test size(F) == (n_z, n_μ, n_μ)
    # Boundary rows and same-direction transfers are unscaled.
    @test all(F[1, :, :] .== 1) && all(F[n_z, :, :] .== 1)
    @test all(F[:, 1:2, 1:2] .== 1) && all(F[:, 3:4, 3:4] .== 1)
    # Down ↔ up transfers differ from 1 on this grid.
    @test all(F[2:(n_z - 1), 3, 1] .!= 1)

    # Any column-stochastic redistribution P then moves exactly the source row's weight.
    P = [0.4 0.1 0.2 0.3; 0.3 0.5 0.1 0.2; 0.2 0.3 0.6 0.1; 0.1 0.1 0.1 0.4]
    @test vec(sum(P; dims = 1)) ≈ ones(n_μ)
    for iz in 2:(n_z - 1), i2 in 1:n_μ
        @test sum(wrow(iz, i1) * F[iz, i1, i2] * P[i1, i2] for i1 in 1:n_μ) ≈ wrow(iz, i2)
    end

    @test_throws "zero-length cell" AURORA.upwind_cell_ratio([0.0, 1.0, 1.0, 2.0], μ)
end
