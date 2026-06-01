@testitem "make_heating_rate_file produces correct output" setup=[FluxTestSetup] begin
    using AURORA
    using NCDatasets

    AURORA.make_heating_rate_file(FluxTestSetup.savedir)

    outfile = joinpath(FluxTestSetup.savedir, "analysis", "heating_rate.nc")
    @test isfile(outfile)

    NCDataset(outfile, "r") do ds
        @test haskey(ds, "heating_rate")
        @test haskey(ds, "altitude")
        @test haskey(ds, "time")
        @test size(ds["heating_rate"]) == (FluxTestSetup.n_z, FluxTestSetup.n_t)
        @test all(Array(ds["heating_rate"]) .>= 0)
    end
end
