@testitem "AltitudeGrid construction" begin
    grid = AltitudeGrid(80, 700)

    @test grid isa AltitudeGrid
    @test grid isa AURORA.AbstractGrid
    @test grid.bottom == 80
    @test grid.top == 700
    @test length(grid.h) == grid.n
    @test length(grid.Δh) == grid.n - 1
    @test all(grid.Δh .> 0) # ensure grid is strictly increasing
    # the grid lands exactly on the requested limits and on the transition altitude
    @test grid.h[1] == 80e3
    @test grid.h[end] == 700e3
    @test 100e3 in grid.h
    # uniform ~150 m spacing below 100 km
    i_transition = findfirst(==(100e3), grid.h)
    @test all(isapprox.(grid.Δh[1:(i_transition - 1)], 150; atol=1))
    # smooth growth above: consecutive steps grow by at most exp(dz_max / growth_scale),
    # except the topmost step which can absorb up to half a step to land exactly on `top`
    ratios = grid.Δh[2:(end - 1)] ./ grid.Δh[1:(end - 2)]
    @test all(ratios .< exp(10 / 60) + 0.01)
    @test grid.Δh[end] < 1.6 * grid.Δh[end - 1]
end

@testitem "AltitudeGrid with bottom above 100 km" begin
    # Used to silently produce a grid starting at 100 km whatever the requested bottom
    grid = AltitudeGrid(250, 500)

    @test grid.h[1] == 250e3
    @test grid.bottom == 250
    @test grid.h[end] == 500e3
    @test all(grid.Δh .> 0)
    @test grid.n < AltitudeGrid(100, 500).n / 2  # F-region grid is much smaller
end

@testitem "AltitudeGrid scale" begin
    grid    = AltitudeGrid(90, 400)
    coarse  = AltitudeGrid(90, 400; scale=2)
    fine    = AltitudeGrid(90, 400; scale=0.5)

    # scale = 1 reproduces the default grid exactly
    @test AltitudeGrid(90, 400; scale=1).h == grid.h

    # point counts scale roughly inversely with the step size
    @test coarse.n < grid.n < fine.n
    @test isapprox(grid.n / coarse.n, 2; rtol=0.1)
    @test isapprox(fine.n / grid.n, 2; rtol=0.1)

    # scaled grids stay valid, with exact endpoints
    for g in (coarse, fine)
        @test g.h[1] == 90e3
        @test g.h[end] == 400e3
        @test all(g.Δh .> 0)
    end
end

@testitem "AltitudeGrid invalid inputs" begin
    @test_throws ArgumentError AltitudeGrid(700, 80)
    @test_throws ArgumentError AltitudeGrid(80, 700; dz0=0)
    @test_throws ArgumentError AltitudeGrid(80, 700; scale=-1)
end

@testitem "suggest_bottom_altitude" begin
    msis_file = joinpath(@__DIR__, "..", "regression", "reference_results",
                         "msis_20051008-2200_70N-19E.txt")

    # 3 keV electrons fully stop near 113 km on this profile; with the default safety
    # margin the suggestion should land a few km below that
    z3 = suggest_bottom_altitude(3000, msis_file)
    @test 100 <= z3 <= 113

    # higher energies penetrate deeper
    @test suggest_bottom_altitude(10_000, msis_file) < z3

    # larger safety margin pushes the boundary further down
    @test suggest_bottom_altitude(3000, msis_file; safety=4) < z3

    # 100 keV electrons penetrate below the bottom of the MSIS file (85 km):
    # warn and fall back to the lowest available altitude
    z100k = @test_logs (:warn, r"penetrate below") suggest_bottom_altitude(100_000, msis_file)
    @test z100k == 85
end

@testitem "suggest_bottom_altitude AuroraModel wrapper" begin
    msis_file = joinpath(@__DIR__, "..", "regression", "reference_results",
                         "msis_20051008-2200_70N-19E.txt")
    iri_file = joinpath(@__DIR__, "..", "regression", "reference_results",
                        "iri_20051008-2200_70N-19E.txt")

    model = AuroraModel([100, 400], 180:-90:0, 3000, msis_file, iri_file)
    # the wrapper reads E_max and the MSIS file from the model
    @test suggest_bottom_altitude(model) ==
          suggest_bottom_altitude(model.energy_grid.E_max, msis_file)
    @test suggest_bottom_altitude(model; safety=4) ==
          suggest_bottom_altitude(model.energy_grid.E_max, msis_file; safety=4)
end

@testitem "EnergyGrid construction" begin
    grid = EnergyGrid(50_000)

    @test grid isa EnergyGrid
    @test grid isa AURORA.AbstractGrid
    @test grid.E_max == 50_000
    @test length(grid.E_edges) == grid.n + 1
    @test length(grid.E_centers) == grid.n
    @test length(grid.ΔE) == grid.n
    @test grid.E_edges[1] > 0
    # TODO: should we allow E_edges to exceed E_max? or stricly enforce E_centers to be below?
    # right now only E_edges[end -1] is guaranteed to be below E_max, but E_edges[end] and
    # E_centers[end] can exceed it
    @test grid.E_edges[end - 1] ≤ 50_000
    @test all(grid.ΔE .> 0)
end

@testitem "EnergyGrid invalid input" begin
    @test_throws Exception EnergyGrid(2e6) # we do not support energies above 1 MeV
end

@testitem "PitchAngleGrid construction" begin
    grid = PitchAngleGrid(180:-10:0)

    @test grid isa PitchAngleGrid
    @test grid isa AURORA.AbstractGrid
    @test grid.n_beams == 18
    @test length(grid.μ_lims) == 19
    @test length(grid.μ_center) == 18
    @test all(-1 .<= grid.μ_lims .<= 1)
    @test all(-1 .<= grid.μ_center .<= 1)
    @test grid.θ_lims[1] == 180
    @test grid.θ_lims[end] == 0
end

@testitem "PitchAngleGrid invalid inputs" begin
    @test_throws ArgumentError PitchAngleGrid(170:-10:0)
    @test_throws ArgumentError PitchAngleGrid(180:-10:10)
    @test_throws ArgumentError PitchAngleGrid(0:10:180)
end
