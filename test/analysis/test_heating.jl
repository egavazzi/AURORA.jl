@testitem "make_heating_rate_file produces correct output" setup=[FluxTestSetup] begin
    using AURORA
    using MAT: matread

    AURORA.make_heating_rate_file(FluxTestSetup.savedir)

    outfile = joinpath(FluxTestSetup.savedir, "heating_rate.mat")
    @test isfile(outfile)

    out = matread(outfile)
    @test haskey(out, "heating_rate")
    @test haskey(out, "h_atm")
    @test haskey(out, "t")

    @test size(out["heating_rate"]) == (FluxTestSetup.n_z, FluxTestSetup.n_t)
    @test all(out["heating_rate"] .>= 0)
end
