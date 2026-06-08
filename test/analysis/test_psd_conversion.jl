@testitem "Conservative F(v_parallel) reduction" begin
    using AURORA

    μ_lims = [-1.0, -0.2, 0.3, 1.0]
    v = [2.0, 5.0]
    Ie = Array{Float64}(undef, 2, 3, 2, 2)

    for i in axes(Ie, 4), it in axes(Ie, 3), j in axes(Ie, 2), iz in axes(Ie, 1)
        Ie[iz, j, it, i] = 10.0 * iz + 3.0 * j + 2.0 * it + i
    end

    result = AURORA.compute_F(Ie, μ_lims, v)

    conserved_density = dropdims(sum(Ie ./ reshape(v, 1, 1, 1, :); dims=(2, 4)), dims=(2, 4))
    reduced_density = dropdims(sum(result.F .* reshape(result.Δvpar, :, 1, 1); dims=1), dims=1)

    @test size(result.F) == (length(result.vpar_centers), size(Ie, 1), size(Ie, 3))
    @test all(result.Δvpar .> 0)
    @test issorted(result.vpar_edges)
    @test result.vpar_edges[1] == -maximum(v)
    @test result.vpar_edges[end] == maximum(v)
    @test any(iszero, result.vpar_edges)
    @test all(isapprox.(diff(result.vpar_edges), first(diff(result.vpar_edges)); rtol=0, atol=1e-12))
    @test conserved_density ≈ reduced_density rtol=1e-12 atol=1e-12

    negative_bins = findall(<(0.0), result.vpar_centers)
    positive_bins = findall(>(0.0), result.vpar_centers)
    @test any(result.F[negative_bins, :, :] .> 0)
    @test any(result.F[positive_bins, :, :] .> 0)
end

@testitem "Custom v_parallel edges preserve density" begin
    using AURORA

    μ_lims = [-1.0, -0.05, 0.05, 1.0]
    v = [3.0, 6.0]
    Ie = fill(4.0, 1, 3, 1, 2)
    vpar_edges = [-6.0, -3.0, -1.0, 0.0, 1.0, 3.0, 6.0]

    result = AURORA.compute_F(Ie, μ_lims, v; vpar_edges=vpar_edges)

    conserved_density = dropdims(sum(Ie ./ reshape(v, 1, 1, 1, :); dims=(2, 4)), dims=(2, 4))
    reduced_density = dropdims(sum(result.F .* reshape(result.Δvpar, :, 1, 1); dims=1), dims=1)

    @test result.vpar_edges == vpar_edges
    @test conserved_density ≈ reduced_density rtol=1e-12 atol=1e-12
    central_bin = findfirst(x -> x == 0.5, result.vpar_centers)
    @test !isnothing(central_bin)
    @test result.F[central_bin, 1, 1] > 0
end

@testitem "make_psd_file writes f and F" begin
    using AURORA
    using NCDatasets

    tmpdir = mktempdir()
    nc_path = joinpath(tmpdir, "simulation_data.nc")

    E_edges   = [10.0, 20.0, 40.0, 80.0]
    E_centers = [15.0, 30.0, 60.0]
    μ_lims    = [-1.0, 0.0, 1.0]
    t_run     = [0.0, 1.0]
    h_atm     = [100e3, 110e3]

    Nz = length(h_atm)
    nμ = length(μ_lims) - 1
    Nt = length(t_run)
    nE = length(E_centers)

    # Ie [Nz, nμ, Nt, nE]
    Ie = zeros(Float32, Nz, nμ, Nt, nE)
    for i_E in 1:nE, i_t in 1:Nt, i_μ in 1:nμ, i_z in 1:Nz
        Ie[i_z, i_μ, i_t, i_E] = Float32(1.0 + i_z + i_μ + i_t + i_E)
    end

    # Write minimal simulation_data.nc
    # Note: unlimited dimension requires indexed writes, not t_v[:] = ...
    NCDataset(nc_path, "c") do ds
        defDim(ds, "altitude",           Nz)
        defDim(ds, "pitch_angle",        nμ)
        defDim(ds, "energy",             nE)
        defDim(ds, "energy_bounds",      nE + 1)
        defDim(ds, "pitch_angle_bounds", nμ + 1)
        defDim(ds, "time",               Inf)

        alt_v = defVar(ds, "altitude",      Float64, ("altitude",))
        alt_v[:] = h_atm
        en_v  = defVar(ds, "energy",        Float64, ("energy",))
        en_v[:] = E_centers
        ee_v  = defVar(ds, "energy_edges",  Float64, ("energy_bounds",))
        ee_v[:] = E_edges
        ml_v  = defVar(ds, "mu_lims",       Float64, ("pitch_angle_bounds",))
        ml_v[:] = μ_lims
        bw_v  = defVar(ds, "beam_weight",   Float64, ("pitch_angle",))
        bw_v[:] = beam_weight(μ_lims)
        defVar(ds, "time", Float64, ("time",))
        defVar(ds, "Ie",   Float32, ("altitude", "pitch_angle", "time", "energy"))
        ds["time"][1:Nt]             = t_run
        ds["Ie"][:, :, 1:Nt, :]      = Ie
    end

    AURORA.make_psd_file(tmpdir; compute = :both)

    outpath = joinpath(tmpdir, "analysis", "psd.nc")
    @test isfile(outpath)

    NCDataset(outpath, "r") do ds
        @test haskey(ds, "f")
        @test haskey(ds, "F")
        @test haskey(ds, "dvpar")
        @test size(ds["f"]) == (Nz, nμ, Nt, nE)
        @test size(ds["F"], 2) == Nz   # (vpar, altitude, time)
        @test size(ds["F"], 3) == Nt

        # Check density conservation
        v = Array(ds["v"])   # [nE]
        Ie_f64 = Float64.(Ie)
        conserved_density = dropdims(sum(Ie_f64 ./ reshape(v, 1, 1, 1, :); dims = (2, 4)), dims = (2, 4))
        dvpar = Array(ds["dvpar"])
        F = Float64.(Array(ds["F"]))
        reduced_density = dropdims(sum(F .* reshape(dvpar, :, 1, 1); dims = 1), dims = 1)
        @test conserved_density ≈ reduced_density rtol = 1e-5 atol = 1e-5
    end
end
