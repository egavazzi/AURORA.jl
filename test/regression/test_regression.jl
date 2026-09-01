@testitem "AURORA steady-state results" begin
    using NCDatasets
    using LazyArtifacts
    altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-30:0;             # (°) angle-limits for the electron beams
    E_max = 500;                   # (eV) upper limit to the energy grid
    B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

    msis_file = joinpath(@__DIR__, "input_data", "msis_20051008-2200_70N-19E.txt")
    iri_file = joinpath(@__DIR__, "input_data", "iri_20051008-2200_70N-19E.txt")

    ## Build the model
    model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)

    ## Define where to save the results
    savedir = mktempdir()

    ## Define input flux
    flux = InputFlux(FlatSpectrum(1e-2; E_min=E_max - 100); beams=1:2)

    ## Run the simulation
    sim = AuroraSimulation(model, flux, savedir; mode=SteadyStateMode())
    initialize!(sim; force_recompute=true, verbose=false) # force recomputation instead of loading from cache to test regressions
    run!(sim)

    ## Analyze the results
    make_volume_excitation_file(savedir)
    make_column_excitation_file(savedir)

    ## Compare the results, allowing a relative difference of 1e-4 (= 0.01%)
    reference_file = joinpath(artifact"reference_results", "SS", "volume_excitation.nc")
    NCDataset(reference_file, "r") do ds_ref
        NCDataset(joinpath(savedir, "analysis", "volume_excitation.nc"), "r") do ds_new
            channel_diffs = Tuple{String, Float64, String}[]
            for name in filter(name -> startswith(name, "Q"), keys(ds_ref))
                values_ref = Array(ds_ref[name])
                values_new = Array(ds_new[name])
                @test all(isapprox.(values_new, values_ref; rtol = 1e-4, atol = 1e-12))

                rel_diff = abs.(values_new .- values_ref) ./
                           max.(abs.(values_new), abs.(values_ref), eps())
                idx = argmax(rel_diff)
                push!(channel_diffs, (String(name), rel_diff[idx], string(idx)))
            end
            println("\nVolume excitation maximum relative differences:")
            println("Channel                  | Maximum relative difference | Index")
            println("-------------------------|-----------------------------|------")
            for (name, diff, idx) in sort(channel_diffs; by=first)
                println(rpad(name, 24), " | ", lpad(string(diff), 27), " | ", idx)
            end
        end
    end

    column_reference = joinpath(artifact"reference_results", "SS", "column_excitation.nc")
    NCDataset(column_reference, "r") do ds_ref
        NCDataset(joinpath(savedir, "analysis", "column_excitation.nc"), "r") do ds_new
            @test Array(ds_new["time"]) == Array(ds_ref["time"])
            channel_diffs = Tuple{String, Float64, String}[]
            for name in ("I_4278", "I_6730", "I_7774", "I_7774_O", "I_7774_O2",
                         "I_8446", "I_8446_O", "I_8446_O2", "I_O1D", "I_O1S")
                values_new = Array(ds_new[name])
                values_ref = Array(ds_ref[name])
                @test size(values_new) == size(values_ref)
                @test all(isapprox.(values_new, values_ref; rtol = 1e-4, atol = 1e-12))

                rel_diff = abs.(values_new .- values_ref) ./
                           max.(abs.(values_new), abs.(values_ref), eps())
                idx = argmax(rel_diff)
                push!(channel_diffs, (name, rel_diff[idx], string(idx)))
            end
            println("\nColumn excitation maximum relative differences:")
            println("Channel                  | Maximum relative difference | Index")
            println("-------------------------|-----------------------------|------")
            for (name, diff, idx) in sort(channel_diffs; by=first)
                println(rpad(name, 24), " | ", lpad(string(diff), 27), " | ", idx)
            end
        end
    end
end


@testitem "AURORA time-dependent results" begin
    using NCDatasets
    using LazyArtifacts
    altitude_lims = [100, 400];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-30:0;             # (°) angle-limits for the electron beams
    E_max = 500;                   # (eV) upper limit to the energy grid
    B_angle_to_zenith = 13;         # (°) angle between the B-field line and the zenith

    msis_file = joinpath(@__DIR__, "input_data", "msis_20051008-2200_70N-19E.txt")
    iri_file = joinpath(@__DIR__, "input_data", "iri_20051008-2200_70N-19E.txt")

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
    initialize!(sim; force_recompute=true, verbose=false) # force recomputation instead of loading from cache to test regressions
    run!(sim)

    ## Analyze the results
    make_volume_excitation_file(savedir)
    make_column_excitation_file(savedir)

    ## Compare the results, allowing a relative difference of 1e-4 (= 0.01%)
    reference_file = joinpath(artifact"reference_results", "TD", "volume_excitation.nc")
    NCDataset(reference_file, "r") do ds_ref
        NCDataset(joinpath(savedir, "analysis", "volume_excitation.nc"), "r") do ds_new
            channel_diffs = Tuple{String, Float64, String}[]
            for name in filter(name -> startswith(name, "Q"), keys(ds_ref))
                values_ref = Array(ds_ref[name])
                values_new = Array(ds_new[name])
                @test all(isapprox.(values_new, values_ref; rtol = 1e-4, atol = 1e-12))

                rel_diff = abs.(values_new .- values_ref) ./
                           max.(abs.(values_new), abs.(values_ref), eps())
                idx = argmax(rel_diff)
                push!(channel_diffs, (String(name), rel_diff[idx], string(idx)))
            end
            println("\nVolume excitation maximum relative differences:")
            println("Channel                  | Maximum relative difference | Index")
            println("-------------------------|-----------------------------|------")
            for (name, diff, idx) in sort(channel_diffs; by=first)
                println(rpad(name, 24), " | ", lpad(string(diff), 27), " | ", idx)
            end
        end
    end

    column_reference = joinpath(artifact"reference_results", "TD", "column_excitation.nc")
    NCDataset(column_reference, "r") do ds_ref
        NCDataset(joinpath(savedir, "analysis", "column_excitation.nc"), "r") do ds_new
            @test Array(ds_new["time"]) == Array(ds_ref["time"])
            channel_diffs = Tuple{String, Float64, String}[]
            for name in ("I_4278", "I_6730", "I_7774", "I_7774_O", "I_7774_O2",
                         "I_8446", "I_8446_O", "I_8446_O2", "I_O1D", "I_O1S")
                values_new = Array(ds_new[name])
                values_ref = Array(ds_ref[name])
                @test size(values_new) == size(values_ref)
                @test all(isapprox.(values_new, values_ref; rtol = 1e-4, atol = 1e-12))

                rel_diff = abs.(values_new .- values_ref) ./
                           max.(abs.(values_new), abs.(values_ref), eps())
                idx = argmax(rel_diff)
                push!(channel_diffs, (name, rel_diff[idx], string(idx)))
            end
            println("\nColumn excitation maximum relative differences:")
            println("Channel                  | Maximum relative difference | Index")
            println("-------------------------|-----------------------------|------")
            for (name, diff, idx) in sort(channel_diffs; by=first)
                println(rpad(name, 24), " | ", lpad(string(diff), 27), " | ", idx)
            end
        end
    end
end
