using MAT

# Here we test that the O1S volume emission rate is the same as some reference results
@testset "AURORA steady-state results" begin
    altitude_lims = [100, 500];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-10:0;             # (°) angle-limits for the electron beams
    E_max = 1000;                   # (eV) upper limit to the energy grid
    B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

    msis_file = find_nrlmsis_file(
        year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
        );
    iri_file = find_iri_file(
        year=2005, month=10, day=8, hour=22, minute=0, lat=70, lon=19, height=85:1:700
        );

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
        msis_file, iri_file, savedir, INPUT_OPTIONS);

    ## Analyze the results
    make_volume_excitation_file(savedir)

    ## Compare the results
    reference_file = "reference_results/Qzt_all_L.mat"
    data_ref = matread(reference_file)
    data_new = matread(joinpath(savedir, "Qzt_all_L.mat"))
    @test isapprox(data_new["QO1S"], data_ref["QO1S"], rtol = 1e-4)

    rm("temp_results", recursive=true)
end
