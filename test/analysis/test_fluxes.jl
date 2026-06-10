# Run one minimal TD simulation once and share the results across all test items.
# θ_lims = 180:-90:0 → beam 1: μ ∈ [-1, 0] (downward), beam 2: μ ∈ [0, 1] (upward).
@testmodule FluxTestSetup begin
    using AURORA
    using NCDatasets

    const savedir = mktempdir()

    let
        altitude_lims = [100, 400]
        θ_lims = 180:-90:0  # 2 beams: beam 1 down (μ<0), beam 2 up (μ>0)
        E_max = 100
        B_angle_to_zenith = 13

        model = AuroraModel(altitude_lims, θ_lims, E_max,
                            find_msis_file(), find_iri_file(),
                            B_angle_to_zenith)
        flux = InputFlux(FlatSpectrum(1.0; E_min=50.0), ConstantModulation();
                         beams=1, z_source=500.0)
        sim = AuroraSimulation(model, flux, savedir;
                               mode=TimeDependentMode(duration = 0.1, dt = 0.01,
                                                      CFL_number = 128, n_loop = 1))
        run!(sim)
    end

    # Read back grid dimensions from simulation_data.nc
    const n_z, n_μ, n_t, n_E = NCDataset(joinpath(savedir, "simulation_data.nc"), "r") do ds
        size(ds["Ie"])   # (altitude, pitch_angle, time, energy)
    end
end

@testitem "make_Ie_top_file extracts top boundary" setup=[FluxTestSetup] begin
    using AURORA
    using NCDatasets

    AURORA.make_Ie_top_file(FluxTestSetup.savedir)

    outfile = joinpath(FluxTestSetup.savedir, "analysis", "Ie_top.nc")
    @test isfile(outfile)

    NCDataset(outfile, "r") do ds
        @test haskey(ds, "Ie_top_raw")
        @test haskey(ds, "Ie_top")

        # Shape: (pitch_angle, time, energy)
        @test size(ds["Ie_top_raw"]) == (FluxTestSetup.n_μ, FluxTestSetup.n_t, FluxTestSetup.n_E)

        # Downward beam (beam 1) received non-zero input flux → positive values at the top
        Ie_top_raw = Array(ds["Ie_top_raw"])
        @test all(Ie_top_raw[1, :, :] .>= 0)
        @test any(Ie_top_raw[1, :, :] .> 0)

        # Upward beam (beam 2) has no direct input — only backscatter, which is non-negative
        @test all(Ie_top_raw[2, :, :] .>= 0)
    end
end

@testitem "make_current_file computes currents" setup=[FluxTestSetup] begin
    using AURORA
    using NCDatasets

    AURORA.make_current_file(FluxTestSetup.savedir)

    outfile = joinpath(FluxTestSetup.savedir, "analysis", "currents.nc")
    @test isfile(outfile)

    NCDataset(outfile, "r") do ds
        @test haskey(ds, "J_up")
        @test haskey(ds, "J_down")
        @test haskey(ds, "IeE_up")
        @test haskey(ds, "IeE_down")
        @test size(ds["J_up"])   == (FluxTestSetup.n_z, FluxTestSetup.n_t)
        @test size(ds["J_down"]) == (FluxTestSetup.n_z, FluxTestSetup.n_t)
        # Downward beam has input → non-zero downward current
        @test all(Array(ds["J_down"]) .>= 0)
        # Backscattered upward flux is also non-negative
        @test all(Array(ds["J_up"]) .>= 0)
    end
end
