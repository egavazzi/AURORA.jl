@testitem "load_volume_excitation" setup=[SharedSimResults] begin
    # Steady-state
    vol_ss = load_volume_excitation(SharedSimResults.ss_dir)
    @test vol_ss isa VolumeExcitationResult
    @test vol_ss.h_atm isa Vector{Float64}
    @test vol_ss.t isa Vector{Float64}
    n_z_ss, n_t_ss = length(vol_ss.h_atm), length(vol_ss.t)
    @test size(vol_ss.Q4278) == (n_z_ss, n_t_ss)
    @test size(vol_ss.QO1S) == (n_z_ss, n_t_ss)
    @test vol_ss.savedir == SharedSimResults.ss_dir

    # Time-dependent
    vol_td = load_volume_excitation(SharedSimResults.td_dir)
    @test vol_td isa VolumeExcitationResult
    n_z_td, n_t_td = length(vol_td.h_atm), length(vol_td.t)
    @test size(vol_td.Q4278) == (n_z_td, n_t_td)
    @test n_t_td > n_t_ss
    @test vol_td.savedir == SharedSimResults.td_dir
end

@testitem "load_column_excitation" setup=[SharedSimResults] begin
    # Steady-state
    col_ss = load_column_excitation(SharedSimResults.ss_dir)
    @test col_ss isa ColumnExcitationResult
    @test col_ss.t isa Vector{Float64}
    @test length(col_ss.I_4278) == length(col_ss.t)
    @test length(col_ss.I_O1S) == length(col_ss.t)

    # Time-dependent
    col_td = load_column_excitation(SharedSimResults.td_dir)
    @test col_td isa ColumnExcitationResult
    @test length(col_td.I_4278) == length(col_td.t)
    @test length(col_td.t) > length(col_ss.t)
end

@testitem "load_results selectors and guard" setup=[SharedSimResults] begin
    full = load_results(SharedSimResults.td_dir)
    @test full isa SimulationResult
    n_z, n_μ, n_t, n_E = size(full.Ie)
    @test length(full.t) == n_t
    @test length(full.h_atm) == n_z
    @test length(full.E_centers) == n_E
    @test length(full.E_edges) == n_E + 1
    @test length(full.μ_lims) == n_μ + 1
    @test length(full.beam_weights) == n_μ

    # Time selector
    sub_t = load_results(SharedSimResults.td_dir; tidx = 1:2)
    @test size(sub_t.Ie, 3) == 2
    @test sub_t.t == full.t[1:2]
    @test sub_t.Ie == full.Ie[:, :, 1:2, :]

    # Energy selector keeps edges / ΔE consistent
    sub_e = load_results(SharedSimResults.td_dir; eidx = 1:3)
    @test size(sub_e.Ie, 4) == 3
    @test sub_e.E_centers == full.E_centers[1:3]
    @test length(sub_e.E_edges) == 4
    @test sub_e.E_edges == full.E_edges[1:4]
    @test sub_e.ΔE == diff(full.E_edges[1:4])

    # Pitch-angle selector keeps μ_lims / beam_weights consistent
    sub_μ = load_results(SharedSimResults.td_dir; μidx = 1:1)
    @test size(sub_μ.Ie, 2) == 1
    @test length(sub_μ.μ_lims) == 2
    @test sub_μ.beam_weights == full.beam_weights[1:1]

    # Guard: a tiny budget throws, Inf bypasses
    @test_throws ArgumentError load_results(SharedSimResults.td_dir; max_bytes = 1)
    big = load_results(SharedSimResults.td_dir; max_bytes = Inf)
    @test size(big.Ie) == size(full.Ie)

    # Non-contiguous selectors are rejected
    @test_throws ArgumentError load_results(SharedSimResults.td_dir; tidx = 1:2:5)
end

@testitem "foreach_Ie_time_chunk streams the full Ie" setup=[SharedSimResults] begin
    full = load_results(SharedSimResults.td_dir)
    n_z, n_μ, n_t, n_E = size(full.Ie)

    # A tiny budget forces multiple chunks.
    reassembled = zeros(n_z, n_μ, n_t, n_E)
    n_chunks = Ref(0)
    foreach_Ie_time_chunk(SharedSimResults.td_dir; max_bytes = 1) do Ie_chunk, t_range
        reassembled[:, :, t_range, :] .= Ie_chunk
        n_chunks[] += 1
    end
    @test n_chunks[] > 1
    @test reassembled == full.Ie

    # An unbounded budget streams everything in a single chunk, same data.
    once = zeros(n_z, n_μ, n_t, n_E)
    n_once = Ref(0)
    foreach_Ie_time_chunk(SharedSimResults.td_dir; max_bytes = Inf) do Ie_chunk, t_range
        once[:, :, t_range, :] .= Ie_chunk
        n_once[] += 1
    end
    @test n_once[] == 1
    @test once == full.Ie
end

@testitem "load_Ie_top" setup=[SharedSimResults] begin
    # Time-dependent
    inp_td = load_Ie_top(SharedSimResults.td_dir)
    @test inp_td isa IeTopResult
    @test inp_td.Ietop isa Array{Float64, 3}
    @test size(inp_td.Ietop, 1) == 2   # 2 beams

    # Steady-state
    inp_ss = load_Ie_top(SharedSimResults.ss_dir)
    @test inp_ss isa IeTopResult
    @test inp_ss.Ietop isa Array{Float64, 3}
    @test size(inp_ss.Ietop, 1) == 2   # 2 beams
    @test size(inp_ss.Ietop, 2) == 1   # 1 time step
end
