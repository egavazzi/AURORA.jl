@testitem "AURORA steady-state results" begin
    using MAT
    altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-30:0;             # (°) angle-limits for the electron beams
    E_max = 500;                   # (eV) upper limit to the energy grid
    B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

    msis_file = joinpath(@__DIR__, "reference_results", "msis_20051008-2200_70N-19E.txt")
    iri_file = joinpath(@__DIR__,  "reference_results", "iri_20051008-2200_70N-19E.txt")

    ## Build the model
    model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)

    ## Define where to save the results
    savedir = mktempdir()

    ## Define input flux
    flux = InputFlux(FlatSpectrum(1e-2; E_min=E_max - 100); beams=1:2)

    ## Run the simulation
    sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode())
    initialize!(sim; force_recompute=true) # force recomputation instead of loading from cache to test regressions
    run!(sim)

    ## Analyze the results
    make_Ie_top_file(savedir)
    make_volume_excitation_file(savedir)
    make_column_excitation_file(savedir)
    make_current_file(savedir)

    ## Compare the results, allowing a relative difference of 1e-4 (= 0.01%)
    reference_file = joinpath(@__DIR__, "reference_results", "SS", "Qzt_all_L.mat")
    data_ref = matread(reference_file)
    data_new = matread(joinpath(savedir, "Qzt_all_L.mat"))
    @test all(isapprox.(data_new["QO1S"], data_ref["QO1S"], rtol = 1e-4))

    ## Print the actual maximum relative difference
    rel_diff = abs.(data_new["QO1S"] .- data_ref["QO1S"]) ./
               max.(abs.(data_new["QO1S"]), abs.(data_ref["QO1S"]), eps())
    println("Maximum relative difference: ", maximum(rel_diff))
end


@testitem "AURORA time-dependent results" begin
    using MAT
    altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-30:0;             # (°) angle-limits for the electron beams
    E_max = 500;                   # (eV) upper limit to the energy grid
    B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

    msis_file = joinpath(@__DIR__, "reference_results", "msis_20051008-2200_70N-19E.txt")
    iri_file = joinpath(@__DIR__, "reference_results", "iri_20051008-2200_70N-19E.txt")

    ## Build the model
    model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)

    ## Define where to save the results
    savedir = mktempdir()

    ## Define input flux
    flux = InputFlux(FlatSpectrum(1e-2; E_min=100.0), SinusoidalFlickering(5.0);
                     beams=1, z_source=1000.0)

    ## Run the simulation
    sim = AuroraSimulation(model, flux, savedir;
                           mode=TimeDependentMode(duration = 0.2, dt = 0.01,
                                                  CFL_number = 128, n_loop = 2))
    initialize!(sim; force_recompute=true) # force recomputation instead of loading from cache to test regressions
    run!(sim)

    ## Analyze the results
    make_Ie_top_file(savedir)
    make_volume_excitation_file(savedir)
    make_column_excitation_file(savedir)
    make_current_file(savedir)

    ## Compare the results, allowing a relative difference of 1e-4 (= 0.01%)
    reference_file = joinpath(@__DIR__, "reference_results", "TD", "Qzt_all_L.mat")
    data_ref = matread(reference_file)
    data_new = matread(joinpath(savedir, "Qzt_all_L.mat"))
    @test all(isapprox.(data_new["QO1S"], data_ref["QO1S"], rtol = 1e-4))

    ## Print the actual maximum relative difference
    rel_diff = abs.(data_new["QO1S"] .- data_ref["QO1S"]) ./
               max.(abs.(data_new["QO1S"]), abs.(data_ref["QO1S"]), eps())
    println("Maximum relative difference: ", maximum(rel_diff))
end
