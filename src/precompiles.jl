using PrecompileTools: @setup_workload, @compile_workload

let
    @setup_workload begin
        energy_grid = EnergyGrid(3000)
        θ_lims = 180:-45:0

        # Atmosphere fixtures shipped with the package, used to compile the full
        # run!/analysis pipeline below. Guarded so precompilation never fails if they
        # are missing.
        fixtures = joinpath(pkgdir(AURORA), "test", "regression", "reference_results")
        msis_file = joinpath(fixtures, "msis_20051008-2200_70N-19E.txt")
        iri_file = joinpath(fixtures, "iri_20051008-2200_70N-19E.txt")
        have_fixtures = isfile(msis_file) && isfile(iri_file)

        @compile_workload begin
            load_cross_sections(energy_grid)
            load_or_compute_scattering(θ_lims, 10; verbose = false)
            cosd.(θ_lims);
            mu_avg(θ_lims);
            beam_weight(θ_lims);

            if have_fixtures
                try
                    model = AuroraModel([100, 150], 180:-90:0, 50, msis_file, iri_file, 13)

                    # Steady-state solve + the analysis files built on top of it
                    ss_dir = mktempdir()
                    ss_flux = InputFlux(FlatSpectrum(1e-2; E_min = 25.0); beams = 1:2)
                    run!(AuroraSimulation(model, ss_flux, ss_dir; mode = SteadyStateMode());
                         verbose = false)
                    make_volume_excitation_file(ss_dir)
                    make_column_excitation_file(ss_dir)
                    make_Ie_top_file(ss_dir)
                    make_current_file(ss_dir)
                    make_heating_rate_file(ss_dir)
                    make_psd_file(ss_dir)

                    td_flux = InputFlux(FlatSpectrum(1e-2; E_min = 25.0),
                                        SinusoidalFlickering(5.0); beams = 1:2)

                    # Multi-step steady state
                    run!(AuroraSimulation(model, td_flux, mktempdir();
                                          mode = SteadyStateMode(duration = 0.02, dt = 0.01));
                         verbose = false)

                    # Crank-Nicolson time-dependent solve
                    run!(AuroraSimulation(model, td_flux, mktempdir();
                                          mode = TimeDependentMode(duration = 0.02, dt = 0.01,
                                                                   CFL_number = 1024,
                                                                   n_loop = 1));
                         verbose = false)
                catch
                    # A failure here must never break precompilation.
                end
            end
        end
    end
end
