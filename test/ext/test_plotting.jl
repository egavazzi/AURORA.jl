@testitem "Plotting extension tests" setup=[SharedSimResults] begin
    using Test
    using AURORA
    using CairoMakie

    CairoMakie.activate!()

    vol = load_volume_excitation(SharedSimResults.td_dir)
    col = load_column_excitation(SharedSimResults.td_dir)
    inp = load_input(SharedSimResults.td_dir)

    @testset "plot_excitation! smoke" begin
        # default call → heatmap
        fig1 = Figure()
        ax1 = Axis(fig1[1, 1])
        @test plot_excitation!(ax1, vol) isa Makie.Heatmap

        # show_contours = true branch
        fig2 = Figure()
        ax2 = Axis(fig2[1, 1])
        @test plot_excitation!(ax2, vol; show_contours = true) isa Makie.Heatmap

        # time_index selects a single profile → Lines
        fig3 = Figure()
        ax3 = Axis(fig3[1, 1]; xscale = log10)
        mid_t = cld(length(vol.t), 2)  # avoid t=0 where Q can be all zeros
        @test plot_excitation!(ax3, vol; time_index = mid_t) isa Makie.Lines

        # out-of-range time_index → ArgumentError from _select_time_index
        fig4 = Figure()
        ax4 = Axis(fig4[1, 1]; xscale = log10)
        @test_throws ArgumentError plot_excitation!(ax4, vol; time_index = 10000)
    end

    @testset "plot_input! smoke" begin
        # scalar beam index
        fig1 = Figure()
        ax1 = Axis(fig1[1, 1]; yscale = log10)
        @test plot_input!(ax1, inp; beams = 1) isa Makie.Heatmap

        # vector beam index
        fig2 = Figure()
        ax2 = Axis(fig2[1, 1]; yscale = log10)
        @test plot_input!(ax2, inp; beams = [1]) isa Makie.Heatmap
    end

    @testset "plot_column_excitation! smoke" begin
        # all wavelengths
        fig1 = Figure()
        ax1 = Axis(fig1[1, 1]; yscale = log10)
        col_plots = plot_column_excitation!(ax1, col)
        @test col_plots isa Vector{Makie.Lines}
        @test length(col_plots) == 6

        # wavelength subset
        fig2 = Figure()
        ax2 = Axis(fig2[1, 1]; yscale = log10)
        result = plot_column_excitation!(ax2, col; wavelengths = [:I_4278])
        @test result isa Vector{Makie.Lines}
        @test length(result) == 1
    end

    @testset "plot_input(sim) smoke" begin
        altitude_lims = [100, 200]
        θ_lims = 180:-45:0
        E_max = 100
        B_angle_to_zenith = 13
        msis_file = find_msis_file()
        iri_file = find_iri_file()
        savedir = "fake_dir"  # won't be used since we're not running the sim

        model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)
        flux = InputFlux(FlatSpectrum(1.0; E_min = 50.0); beams = 1:2)
        sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode())

        @test plot_input(sim) isa Figure
    end

    @testset "multi-step SS smoke" begin
        vol_ms = load_volume_excitation(SharedSimResults.ms_ss_dir)
        col_ms = load_column_excitation(SharedSimResults.ms_ss_dir)
        inp_ms = load_input(SharedSimResults.ms_ss_dir)

        # plot_excitation! heatmap over time
        fig1 = Figure()
        ax1 = Axis(fig1[1, 1])
        @test plot_excitation!(ax1, vol_ms) isa Makie.Heatmap

        # plot_excitation! single time slice
        fig2 = Figure()
        ax2 = Axis(fig2[1, 1]; xscale = log10)
        mid_t = cld(length(vol_ms.t), 2)
        @test plot_excitation!(ax2, vol_ms; time_index = mid_t) isa Makie.Lines

        # plot_column_excitation!
        fig3 = Figure()
        ax3 = Axis(fig3[1, 1]; yscale = log10)
        @test plot_column_excitation!(ax3, col_ms) isa Vector{Makie.Lines}

        # plot_input! from loaded data
        fig4 = Figure()
        ax4 = Axis(fig4[1, 1]; yscale = log10)
        @test plot_input!(ax4, inp_ms; beams = 1) isa Makie.Heatmap

        # plot_input from sim object
        altitude_lims = [100, 200]
        θ_lims = 180:-45:0
        E_max = 100
        B_angle_to_zenith = 13
        msis_file = find_msis_file()
        iri_file = find_iri_file()
        savedir = "fake_dir"  # won't be used since we're not running the sim

        model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)
        flux = InputFlux(FlatSpectrum(1.0; E_min = 50.0), SinusoidalFlickering(5.0); beams = 1:2)
        sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode(0.04, 0.01))
        @test plot_input(sim) isa Figure
    end

    # This testset will display figures when run on a personal machine, which is to be
    # expected as we have `display(fig)` internally. We just have to live with it.
    @testset "animation smoke" begin
        # test the plot_input=false branch
        fig = animate_Ie_in_time(SharedSimResults.td_dir;
                                 save_to_file = false,
                                 plot_input = false,
                                 dt_steps = 5)
        @test fig isa Figure

        # test the plot_input=true branch
        fig = animate_Ie_in_time(SharedSimResults.td_dir;
                                 save_to_file = false,
                                 plot_input = true,
                                 input_angle_cone = [90, 180],
                                 dt_steps = 5)
        @test fig isa Figure

        # dt_steps < 1 must throw
        @test_throws Exception animate_Ie_in_time(SharedSimResults.td_dir;
            save_to_file = false, dt_steps = 0)
    end

    @testset "plot_model smoke" begin
        model = SharedSimResults.model

        # default: all panels
        figs = plot_model(model)
        @test figs isa Dict{Symbol, Figure}
        @test Set(keys(figs)) == Set([:atmosphere, :energy_levels, :energy_grid,
                                      :cross_sections, :phase_functions, :scattering])
        for (_, fig) in figs
            @test fig isa Figure
        end

        # panels kwarg: subset
        figs2 = plot_model(model; panels = [:atmosphere, :energy_grid])
        @test Set(keys(figs2)) == Set([:atmosphere, :energy_grid])

        # single panel
        figs3 = plot_model(model; panels = [:cross_sections])
        @test haskey(figs3, :cross_sections)
        @test length(figs3) == 1

        # unknown panel warns but doesn't error
        figs4 = @test_logs (:warn, r"Unknown panel") plot_model(model; panels = [:bogus])
        @test isempty(figs4)
    end

end
