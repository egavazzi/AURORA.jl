@testitem "AURORA steady-state results" begin
    using MAT
    altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-30:0;             # (°) angle-limits for the electron beams
    E_max = 500;                   # (eV) upper limit to the energy grid
    B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

    msis_file = "reference_results/msis_20051008-2200_70N-19E.txt"
    iri_file = "reference_results/iri_20051008-2200_70N-19E.txt"

    ## Define where to save the results
    root_savedir = "temp_results/"   # name of the root folder
    name_savedir = "temp/"   # name of the experiment folder
    savedir = make_savedir(root_savedir, name_savedir; behavior = "custom")

    ## Define input parameters
    input_type = "constant_onset"
    IeE_tot = 1e-2;             # (W/m²) total energy flux of the FAB
    z₀ = altitude_lims[2];      # (km) altitude of the source
    E_min = E_max - 100;        # (eV) bottom energy of the FAB
    Beams = 1:2;                # beam numbers for the precipitation, starting with field aligned down
    t0 = 0;                     # (s) time of start for smooth transition
    t1 = 0;                     # (s) time of end for smooth transition
    INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, Beams, t0, t1);

    ## Run the simulation
    calculate_e_transport_steady_state(altitude_lims, θ_lims, E_max, B_angle_to_zenith,
                                       msis_file, iri_file, savedir, INPUT_OPTIONS)

    ## Analyze the results
    make_Ie_top_file(savedir)
    make_volume_excitation_file(savedir)
    make_column_excitation_file(savedir)
    make_current_file(savedir)

    ## Compare the results, allowing a relative difference of 1e-4 (= 0.01%)
    reference_file = "reference_results/SS/Qzt_all_L.mat"
    data_ref = matread(reference_file)
    data_new = matread(joinpath(savedir, "Qzt_all_L.mat"))
    @test all(isapprox.(data_new["QO1S"], data_ref["QO1S"], rtol = 1e-4))

    ## Print the actual maximum relative difference
    rel_diff = abs.(data_new["QO1S"] .- data_ref["QO1S"]) ./
               max.(abs.(data_new["QO1S"]), abs.(data_ref["QO1S"]), eps())
    println("Maximum relative difference: ", maximum(rel_diff))

    rm("temp_results", recursive=true)
end


@testitem "AURORA time-dependent results" begin
    using MAT
    altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-30:0;             # (°) angle-limits for the electron beams
    E_max = 500;                   # (eV) upper limit to the energy grid
    B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

    t_total = 0.2
    dt = 0.01
    CFL_number = 128;

    msis_file = "reference_results/msis_20051008-2200_70N-19E.txt"
    iri_file = "reference_results/iri_20051008-2200_70N-19E.txt"

    ## Define where to save the results
    root_savedir = "temp_results/"   # name of the root folder
    name_savedir = "temp/"   # name of the experiment folder
    savedir = make_savedir(root_savedir, name_savedir; behavior = "custom")

    ## Define input parameters
    input_type = "flickering";
    IeE_tot = 1e-2;             # (W/m²) total energy flux of the FAB
    z₀ = 1000;                  # (km) altitude of the source
    E_min = 100;                # (eV) bottom energy of the FAB
    f = 5;                      # (Hz) frequence of the modulation
    Beams = 1;                  # beam numbers for the precipitation, starting with field aligned down
    modulation = "sinus";       # type of the modulation ("square" or "sinus")
    INPUT_OPTIONS = (;input_type, IeE_tot, z₀, E_min, f, Beams, modulation);

    ## Run the simulation
    calculate_e_transport(altitude_lims, θ_lims, E_max, B_angle_to_zenith, t_total, dt,
                          msis_file, iri_file, savedir, INPUT_OPTIONS, CFL_number; n_loop = 2)

    ## Analyze the results
    make_Ie_top_file(savedir)
    make_volume_excitation_file(savedir)
    make_column_excitation_file(savedir)
    make_current_file(savedir)

    ## Compare the results, allowing a relative difference of 1e-4 (= 0.01%)
    reference_file = "reference_results/TD/Qzt_all_L.mat"
    data_ref = matread(reference_file)
    data_new = matread(joinpath(savedir, "Qzt_all_L.mat"))
    @test all(isapprox.(data_new["QO1S"], data_ref["QO1S"], rtol = 1e-4))

    ## Print the actual maximum relative difference
    rel_diff = abs.(data_new["QO1S"] .- data_ref["QO1S"]) ./
               max.(abs.(data_new["QO1S"]), abs.(data_ref["QO1S"]), eps())
    println("Maximum relative difference: ", maximum(rel_diff))

    rm("temp_results", recursive=true)
end
