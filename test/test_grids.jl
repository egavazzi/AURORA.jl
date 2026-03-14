@testitem "AltitudeGrid construction" begin
    grid = AltitudeGrid(80, 700)

    @test grid isa AltitudeGrid
    @test grid isa AURORA.AbstractGrid
    @test grid.bottom == 80
    @test grid.top == 700
    @test length(grid.h) == grid.n
    @test length(grid.Δh) == grid.n - 1
    @test all(grid.Δh .> 0) # ensure grid is strictly increasing
    @test grid.h[1] ≈ 80e3
    @test grid.h[end] ≤ 700e3
end

@testitem "AltitudeGrid invalid inputs" begin
    @test_throws Exception AltitudeGrid(700, 80)
end

@testitem "EnergyGrid construction" begin
    grid = EnergyGrid(50_000)

    @test grid isa EnergyGrid
    @test grid isa AURORA.AbstractGrid
    @test grid.E_max == 50_000
    @test length(grid.E) == grid.n
    @test length(grid.dE) == grid.n
    @test grid.E[1] > 0
    @test grid.E[end] ≤ 50_000
    @test all(grid.dE .> 0) # ensure grid is strictly increasing
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
