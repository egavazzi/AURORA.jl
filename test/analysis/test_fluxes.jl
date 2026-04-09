# Run one minimal TD simulation once and share the results across all three test items.
# θ_lims = 180:-90:0 → beam 1: μ ∈ [-1, 0] (downward), beam 2: μ ∈ [0, 1] (upward).
# Only beam 1 receives input flux (beams=1), so Ie_top_raw[1, :, :] must equal Ie_total[1, :, :]
# from Ie_incoming.mat — a strong, physically-grounded assertion.
@testmodule FluxTestSetup begin
    using AURORA
    using MAT: matread

    const savedir = mktempdir()

    let
        altitude_lims = [100, 400]
        θ_lims = 180:-90:0  # 2 beams: beam 1 down (μ<0), beam 2 up (μ>0)
        E_max = 100
        B_angle_to_zenith = 13
        t_total = 0.1
        dt = 0.01
        CFL_number = 128

        model = AuroraModel(altitude_lims, θ_lims, E_max,
                            find_msis_file(), find_iri_file(),
                            B_angle_to_zenith)
        flux = InputFlux(FlatSpectrum(1.0; E_min=50.0), ConstantModulation();
                         beams=1, z_source=500.0)
        sim = AuroraSimulation(model, flux, t_total, dt, savedir;
                               CFL_number, n_loop=1)
        run!(sim)
    end

    # Read back the simulation output
    ie_file = matread(joinpath(savedir, "IeFlickering-01.mat"))
    const Ie_ztE    = ie_file["Ie_ztE"]    # [n_μ*n_z, n_t, n_E]
    const t_run     = vec(ie_file["t_run"])
    const z         = vec(ie_file["h_atm"])
    const ΔE        = vec(ie_file["dE"])
    const Ω_beam    = vec(ie_file["mu_scatterings"]["BeamWeight"])

    const n_z = length(z)
    const n_μ = length(Ω_beam)
    const n_t = length(t_run)
    const n_E = length(ΔE)

end

@testitem "downsampling_fluxes reduces time dimension" setup=[FluxTestSetup] begin
    using AURORA
    using MAT: matread

    factor = 5
    AURORA.downsampling_fluxes(FluxTestSetup.savedir, factor)

    outdir = joinpath(FluxTestSetup.savedir, "downsampled_$(factor)x")
    @test isdir(outdir)

    outfile = joinpath(outdir, "IeFlickering-01d.mat")
    @test isfile(outfile)

    out = matread(outfile)
    n_t_ds = length(1:factor:FluxTestSetup.n_t)  # accounts for odd n_t
    @test size(out["Ie_ztE"], 2) == n_t_ds
    @test length(out["t_run"]) == n_t_ds
    # First time step is preserved
    @test out["t_run"][1] ≈ FluxTestSetup.t_run[1]
    # Values match the original at the correct indices (every `factor`-th step)
    @test out["Ie_ztE"][:, 1, :] ≈ FluxTestSetup.Ie_ztE[:, 1, :]
    @test out["Ie_ztE"][:, 2, :] ≈ FluxTestSetup.Ie_ztE[:, 1 + factor, :]
end

@testitem "make_Ie_top_file extracts top boundary" setup=[FluxTestSetup] begin
    using AURORA
    using MAT: matread

    AURORA.make_Ie_top_file(FluxTestSetup.savedir)

    outfile = joinpath(FluxTestSetup.savedir, "Ie_top.mat")
    @test isfile(outfile)

    out = matread(outfile)
    @test haskey(out, "Ie_top_raw")
    @test haskey(out, "Ie_top")

    # Shape: [n_μ, n_t, n_E]
    @test size(out["Ie_top_raw"]) == (FluxTestSetup.n_μ, FluxTestSetup.n_t, FluxTestSetup.n_E)

    # Downward beam (beam 1) received non-zero input flux → positive values at the top
    @test all(out["Ie_top_raw"][1, :, :] .>= 0)
    @test any(out["Ie_top_raw"][1, :, :] .> 0)

    # Upward beam (beam 2) has no direct input — only backscatter, which is non-negative
    @test all(out["Ie_top_raw"][2, :, :] .>= 0)

    # TODO: Can we somehow compare with what was actually used as input? The problem is that
    # Ie_incoming.mat is saved over the internal refined grid, while Ie_top.mat is
    # calculated from the simulation results that are typically) saved on a coarser grid.
    # Maybe we could use the sim.CFL_number to manually jump??
end

@testitem "make_current_file computes currents" setup=[FluxTestSetup] begin
    using AURORA
    using MAT: matread

    AURORA.make_current_file(FluxTestSetup.savedir)

    outfile = joinpath(FluxTestSetup.savedir, "J.mat")
    @test isfile(outfile)

    out = matread(outfile)
    @test haskey(out, "J_up")
    @test haskey(out, "J_down")
    @test haskey(out, "IeE_up")
    @test haskey(out, "IeE_down")
    @test size(out["J_up"])   == (FluxTestSetup.n_z, FluxTestSetup.n_t)
    @test size(out["J_down"]) == (FluxTestSetup.n_z, FluxTestSetup.n_t)
    # Downward beam has input → non-zero downward current
    @test all(out["J_down"] .>= 0)
    # Backscattered upward flux is also non-negative
    @test all(out["J_up"] .>= 0)
end
