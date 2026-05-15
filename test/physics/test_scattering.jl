@testitem "ScatteringData construction" begin
    grid = PitchAngleGrid(180:-30:0)
    sd = AURORA.ScatteringData(grid)

    @test sd isa AURORA.ScatteringData
    @test size(sd.Ω_subbeam_relative, 1) == grid.n_beams
    @test size(sd.P_scatter, 1) == size(sd.P_scatter, 2)
    @test length(sd.Ω_beam) == grid.n_beams
    @test length(sd.θ_scatter) > 0
end

@testitem "Scattering cache lifecycle" begin
    using JLD2: jldopen

    cache_files(dir) = isdir(dir) ? filter(name -> endswith(name, ".jld2"), readdir(dir)) : String[]
    compatible_cache_count(dir) = count(cache_files(dir)) do filename
        filepath = joinpath(dir, filename)
        jldopen(filepath, "r") do file
            string(file["version_AURORA"]) == AURORA.cache_version_string()
        end
    end

    cache_root = mktempdir()
    θ_lims = 180:-90:0
    n_direction = 10

    save_policy = AURORA.CachePolicy(force_recompute=true, save_cache=true, cache_root=cache_root)
    load_policy = AURORA.CachePolicy(cache_root=cache_root)
    skip_save_policy = AURORA.CachePolicy(force_recompute=true, save_cache=false,
                                          cache_root=joinpath(cache_root, "skip_save"))

    P_scatter, Ω_subbeam_relative, θ_scatter =
        AURORA.load_or_compute_scattering_cache(θ_lims, n_direction; verbose=false, policy=save_policy)

    scattering_dir = joinpath(cache_root, "e_scattering")
    @test length(cache_files(scattering_dir)) == 1

    P_loaded, Ω_loaded, θ_loaded =
        AURORA.load_or_compute_scattering_cache(θ_lims, n_direction; verbose=false, policy=load_policy)
    @test P_loaded == P_scatter
    @test Ω_loaded == Ω_subbeam_relative
    @test θ_loaded == θ_scatter

    scattering_file = joinpath(scattering_dir, only(cache_files(scattering_dir)))
    payload = jldopen(scattering_file, "r") do file
        (
            P_scatter = file["P_scatter"],
            Ω_subbeam_relative = file["Ω_subbeam_relative"],
            theta_scatter = file["theta_scatter"],
            theta_lims = file["theta_lims"],
            n_direction = file["n_direction"],
        )
    end
    rm(scattering_file; force=true)
    jldopen(scattering_file, "w") do file
        file["version_AURORA"] = "0.0.0"
        file["P_scatter"] = payload.P_scatter
        file["Ω_subbeam_relative"] = payload.Ω_subbeam_relative
        file["theta_scatter"] = payload.theta_scatter
        file["theta_lims"] = payload.theta_lims
        file["n_direction"] = payload.n_direction
    end

    AURORA.load_or_compute_scattering_cache(θ_lims, n_direction; verbose=false, policy=load_policy)
    @test compatible_cache_count(scattering_dir) >= 1

    AURORA.load_or_compute_scattering_cache(θ_lims, n_direction; verbose=false, policy=skip_save_policy)
    skip_dir = joinpath(cache_root, "skip_save", "e_scattering")
    @test isempty(cache_files(skip_dir))

    AURORA.clear_scattering_cache!(cache_root=cache_root)
    @test isempty(cache_files(scattering_dir))
end

@testitem "Rotation matrices" begin
    # Define some parameters for the tests
    θ_lims = 180:-10:0
    n_dirs = 720
    P_scatter, Ω_subbeam_relative, θ₁ = AURORA.calculate_scattering_matrices(θ_lims, n_dirs);
    BeamW = beam_weight(θ_lims);

    # check if θ_lims is symmetric
    check = false
    if iseven(length(θ_lims) - 1)
        n_θhalf = Int((length(θ_lims) - 1) / 2)
        check = all(θ_lims[1:n_θhalf + 1] .+ θ_lims[end:-1:n_θhalf + 1] .== 180)
    end

    @testset "Ω_subbeam_relative" begin
        # check normalization
        @test all(sum(Ω_subbeam_relative, dims=2) .≈ 1)
    end

    @testset "P_scatter" begin
        # check normalization (for any i and j, sum(P_scatter[i, j, :], dims=3) = 1)
        @test all(sum(P_scatter, dims = 3) .≈ 1)
        if check
            # if θ_lims is symmetric, check that sum(P_scatter, dims=(1, 2)) is symmetric
            for i in 1:n_θhalf
                @test isapprox(sum(P_scatter[:, :, i]), sum(P_scatter[:, :, end + 1 - i]))
            end
        end
    end
end

@testitem "B2B matrices" begin
    # Define some parameters for the tests
    θ_lims = 180:-10:0
    n_dirs = 720
    P_scatter, Ω_subbeam_relative, θ₁ = AURORA.calculate_scattering_matrices(θ_lims, n_dirs);

    # check if θ_lims is symmetric
    check = false
    if iseven(length(θ_lims) - 1)
        n_θhalf = Int((length(θ_lims) - 1) / 2)
        check = all(θ_lims[1:n_θhalf + 1] .+ θ_lims[end:-1:n_θhalf + 1] .== 180)
    end

    # if θ_lims is symmetric, the B2B matrices should be symmetric
    if check
        _, E_centers, _ = AURORA.make_energy_grid(7000)
        phaseN2e, phaseN2i = AURORA.phase_fcn_N2(θ₁, E_centers);
        phaseO2e, phaseO2i = AURORA.phase_fcn_O2(θ₁, E_centers);
        phaseOe, phaseOi = AURORA.phase_fcn_O(θ₁, E_centers);
        phase_fcn = ((phaseN2e, phaseN2i), (phaseO2e, phaseO2i), (phaseOe, phaseOi));

        B2B_fragment = AURORA.prepare_beams2beams(Ω_subbeam_relative, P_scatter);

        iE = 400 # ~3730 eV (just to pick one)
        for i in 1:3
            phase_fcn_e = AURORA.convert_phase_fcn_to_3D(phase_fcn[i][1][:, iE], θ₁);
            phase_fcn_i = AURORA.convert_phase_fcn_to_3D(phase_fcn[i][2][:, iE], θ₁);
            B2B_elastic = AURORA.beams2beams(phase_fcn_e, B2B_fragment);
            B2B_inelastic = AURORA.beams2beams(phase_fcn_i, B2B_fragment);

            @test all(sum(B2B_elastic, dims = 1) .≈ 1)
            @test all(sum(B2B_inelastic, dims = 1) .≈ 1)
        end
    end
end
