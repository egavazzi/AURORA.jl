using AURORA
using TestItemRunner
using TestItems

@testmodule SharedSimResults begin
    using AURORA

    # Minimal grid shared by both simulations
    altitude_lims = [100, 200]
    θ_lims = 180:-90:0          # 2 beams
    E_max = 100
    B_angle_to_zenith = 13
    msis_file = find_msis_file()
    iri_file = find_iri_file()
    model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)

    # --- Steady-state simulation ---
    ss_dir = mktempdir()
    let flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0); beams = 1:2)
        sim = AuroraSimulation(model, flux, ss_dir; mode=SteadyStateMode())
        run!(sim)
    end
    make_volume_excitation_file(ss_dir)
    make_column_excitation_file(ss_dir)

    # --- Multi-step steady-state simulation ---
    ms_ss_dir = mktempdir()
    let flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0), SinusoidalFlickering(5.0); beams = 1:2)
        sim = AuroraSimulation(model, flux, ms_ss_dir; mode=SteadyStateMode(duration=0.1, dt=0.01))
        run!(sim)
    end
    make_volume_excitation_file(ms_ss_dir)
    make_column_excitation_file(ms_ss_dir)

    # --- Time-dependent simulation (2 loops) ---
    td_dir = mktempdir()
    let flux = InputFlux(FlatSpectrum(1e-2; E_min = 50.0), SinusoidalFlickering(5.0); beams = 1:2)
        sim = AuroraSimulation(model, flux, td_dir;
                               mode=TimeDependentMode(duration=0.1, dt=0.01,
                                                      CFL_number=128, n_loop=2))
        run!(sim)
    end
    make_volume_excitation_file(td_dir)
    make_column_excitation_file(td_dir)
end
@run_package_tests verbose=true
