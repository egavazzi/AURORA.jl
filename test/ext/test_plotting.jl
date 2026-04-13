@testitem "Plotting extension smoke tests" begin
    using Test
    using AURORA
    using CairoMakie
    using MAT: matwrite

    CairoMakie.activate!()

    positive_matrix(n1, n2; scale = 1.0) = reshape(scale .* collect(1.0:(n1 * n2)), n1, n2)
    positive_vector(n; offset = 1.0) = collect(offset:(offset + n - 1))

    function make_volume_result(; n_z = 8, n_t = 4)
        h_atm = collect(range(100e3, 300e3; length = n_z))
        t = collect(range(0.0, 0.3; length = n_t))

        return VolumeExcitationResult(
            positive_matrix(n_z, n_t; scale = 1.0),
            positive_matrix(n_z, n_t; scale = 1.1),
            positive_matrix(n_z, n_t; scale = 1.2),
            positive_matrix(n_z, n_t; scale = 1.3),
            positive_matrix(n_z, n_t; scale = 1.4),
            positive_matrix(n_z, n_t; scale = 1.5),
            positive_matrix(n_z, n_t; scale = 1.6),
            positive_matrix(n_z, n_t; scale = 1.7),
            positive_matrix(n_z, n_t; scale = 1.8),
            positive_matrix(n_z, n_t; scale = 1.9),
            positive_matrix(n_z, n_t; scale = 2.0),
            positive_matrix(n_z, n_t; scale = 2.1),
            positive_matrix(n_z, n_t; scale = 2.2),
            h_atm,
            t,
            nothing,
        )
    end

    function make_column_result(; n_t = 4)
        t = collect(range(0.0, 0.3; length = n_t))

        return ColumnExcitationResult(
            positive_vector(n_t; offset = 1.0),
            positive_vector(n_t; offset = 2.0),
            positive_vector(n_t; offset = 3.0),
            positive_vector(n_t; offset = 4.0),
            positive_vector(n_t; offset = 5.0),
            positive_vector(n_t; offset = 6.0),
            positive_vector(n_t; offset = 7.0),
            positive_vector(n_t; offset = 8.0),
            positive_vector(n_t; offset = 9.0),
            positive_vector(n_t; offset = 10.0),
            t,
        )
    end

    function make_input_result(; n_beams = 2, n_t = 4, n_E = 6)
        Ietop = Array{Float64}(undef, n_beams, n_t, n_E)
        for i_beam in 1:n_beams, i_t in 1:n_t, i_E in 1:n_E
            Ietop[i_beam, i_t, i_E] = i_beam + i_t + i_E
        end

        return IeTopResult(
            Ietop,
            collect(range(0.0, 0.3; length = n_t)),
            [50.0, 100.0, 150.0, 250.0, 400.0, 700.0, 1000.0],
            [75.0, 125.0, 200.0, 325.0, 550.0, 850.0],
            [50.0, 50.0, 100.0, 150.0, 300.0, 300.0],
            [-1.0, 0.0, 1.0],
        )
    end

    function make_animation_data(tmpdir)
        μ_lims = [-1.0, 0.0, 1.0]
        z = [100e3, 150e3, 200e3]
        t_run = collect(range(0.0, 0.1; length = 10))
        E_centers = [100.0, 200.0, 400.0]
        dE = [100.0, 100.0, 200.0]
        n_μ = length(μ_lims) - 1
        n_z = length(z)
        n_t = length(t_run)
        n_E = length(E_centers)

        Ie_ztE = Array{Float64}(undef, n_μ * n_z, n_t, n_E)
        for i in eachindex(Ie_ztE)
            Ie_ztE[i] = 1.0 + i / 10
        end

        matwrite(joinpath(tmpdir, "IeFlickering-1.mat"), Dict(
            "Ie_ztE" => Ie_ztE,
            "mu_lims" => μ_lims,
            "t_run" => t_run,
            "E_centers" => E_centers,
            "h_atm" => z,
            "dE" => dE,
        ))
    end

    @testset "building-block smoke and types" begin
        vol = make_volume_result()
        col = make_column_result()
        input = make_input_result()

        fig_heatmap = Figure()
        ax_heatmap = Axis(fig_heatmap[1, 1])
        @test plot_excitation!(ax_heatmap, vol) isa Makie.Heatmap

        fig_profile = Figure()
        ax_profile = Axis(fig_profile[1, 1]; xscale = log10)
        @test plot_excitation!(ax_profile, vol; time_index = 1) isa Makie.Lines

        fig_input_scalar = Figure()
        ax_input_scalar = Axis(fig_input_scalar[1, 1]; yscale = log10)
        @test plot_input!(ax_input_scalar, input; beams = 1) isa Makie.Heatmap

        fig_input_vector = Figure()
        ax_input_vector = Axis(fig_input_vector[1, 1]; yscale = log10)
        @test plot_input!(ax_input_vector, input; beams = [1]) isa Makie.Heatmap

        fig_column = Figure()
        ax_column = Axis(fig_column[1, 1]; yscale = log10)
        col_plots = plot_column_excitation!(ax_column, col)
        @test col_plots isa Vector{Makie.Lines}
        @test length(col_plots) == 6
    end

    @testset "plot_input(sim) wrapper" begin
        mktempdir() do savedir
            altitude_lims = [100, 200]
            θ_lims = 180:-45:0
            E_max = 100
            B_angle_to_zenith = 13
            msis_file = find_msis_file()
            iri_file = find_iri_file()

            model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)
            flux = InputFlux(FlatSpectrum(1.0; E_min = 50.0); beams = 1:2)
            sim = AuroraSimulation(model, flux, savedir)

            @test plot_input(sim) isa Figure
        end
    end

    # This tests will display a figure, which is to be expected as we have a `display(fig)`
    # internally. We just have to live with it.
    @testset "animation smoke" begin
        mktempdir() do tmpdir
            make_animation_data(tmpdir)

            fig = animate_Ie_in_time(tmpdir; save_to_file = false, plot_input = false, dt_steps = 2)
            @test fig isa Figure
        end
    end
end
