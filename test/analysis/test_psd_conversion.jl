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
    using MAT: matopen, matread

    tmpdir = mktempdir()
    filepath = joinpath(tmpdir, "IeFlickering-01.mat")

    E = [10.0, 20.0, 40.0]
    μ_lims = [-1.0, 0.0, 1.0]
    t_run = [0.0, 1.0]
    h_atm = [100e3, 110e3]

    Nz = length(h_atm)
    nμ = length(μ_lims) - 1
    Nt = length(t_run)
    nE = length(E)

    Ie = zeros(Float64, nμ * Nz, Nt, nE)
    for i_E in 1:nE, i_t in 1:Nt, i_μ in 1:nμ, i_z in 1:Nz
        Ie[(i_μ - 1) * Nz + i_z, i_t, i_E] = 1.0 + i_z + i_μ + i_t + i_E
    end

    matopen(filepath, "w") do io
        write(io, "Ie_ztE", Ie)
        write(io, "E", E)
        write(io, "mu_lims", μ_lims)
        write(io, "t_run", t_run)
        write(io, "h_atm", h_atm)
    end

    AURORA.make_psd_file(tmpdir; compute = :both)

    outpath = joinpath(tmpdir, "psd", "psd-01.mat")
    @test isfile(outpath)

    out = matread(outpath)
    @test haskey(out, "f")
    @test haskey(out, "F")
    @test haskey(out, "dvpar")
    @test size(out["f"]) == (Nz, nμ, Nt, nE)
    @test size(out["F"], 2) == Nz
    @test size(out["F"], 3) == Nt

    Ie_reshaped = reshape(Ie, Nz, nμ, Nt, nE)
    v = out["v"]
    conserved_density = dropdims(sum(Ie_reshaped ./ reshape(v, 1, 1, 1, :); dims = (2, 4)), dims = (2, 4))
    reduced_density = dropdims(sum(out["F"] .* reshape(out["dvpar"], :, 1, 1); dims = 1), dims = 1)
    @test conserved_density ≈ reduced_density rtol = 1e-10 atol = 1e-10
end
