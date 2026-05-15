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
