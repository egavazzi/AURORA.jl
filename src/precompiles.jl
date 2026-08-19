using PrecompileTools: @setup_workload, @compile_workload

let
    @setup_workload begin
        energy_grid = EnergyGrid(3000)
        θ_lims = 180:-45:0

        msis_file = joinpath(pkgdir(AURORA), "test", "regression", "input_data",
                             "msis_20051008-2200_70N-19E.txt")
        iri_file = joinpath(pkgdir(AURORA), "test", "regression", "input_data",
                            "iri_20051008-2200_70N-19E.txt")

        @compile_workload begin
            load_cross_sections(energy_grid)
            load_or_compute_scattering(θ_lims, 10; verbose = false)
            cosd.(θ_lims);
            mu_avg(θ_lims);
            beam_weight(θ_lims);

            model = AuroraModel([100, 150], 180:-90:0, 50, msis_file, iri_file, 13)

            # Steady-state solve + the analysis files built on top of it
            ss_dir = mktempdir()
            ss_flux = InputFlux(FlatSpectrum(1e-2; E_min = 25.0); beams = 1:2)
            mode = SteadyStateMode()
            run!(AuroraSimulation(model, ss_flux, ss_dir; mode), verbose = false)
            make_volume_excitation_file(ss_dir)
            make_column_excitation_file(ss_dir)
            make_Ie_top_file(ss_dir)
            make_current_file(ss_dir)
            make_heating_rate_file(ss_dir)
            make_psd_file(ss_dir)

            # Multi-step steady state
            td_flux = InputFlux(FlatSpectrum(1e-2; E_min = 25.0),
                                SinusoidalFlickering(5.0); beams = 1:2)
            mode = SteadyStateMode(duration = 0.02, dt = 0.01)
            run!(AuroraSimulation(model, td_flux, mktempdir(); mode), verbose = false)

            # Crank-Nicolson time-dependent solve
            td_flux = InputFlux(FlatSpectrum(1e-2; E_min = 25.0),
                                SinusoidalFlickering(5.0); beams = 1:2)
            mode = TimeDependentMode(duration = 0.02, dt = 0.01, CFL_number = 1000,
                                     n_loop = 1)
            run!(AuroraSimulation(model, td_flux, mktempdir(); mode), verbose = false)
        end
    end
end
