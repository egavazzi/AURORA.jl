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

@testitem "load_input" setup=[SharedSimResults] begin
    # Time-dependent
    inp_td = load_input(SharedSimResults.td_dir)
    @test inp_td isa IeTopResult
    @test inp_td.Ietop isa Array{Float64, 3}
    @test size(inp_td.Ietop, 1) == 2   # 2 beams

    # Steady-state
    inp_ss = load_input(SharedSimResults.ss_dir)
    @test inp_ss isa IeTopResult
    @test inp_ss.Ietop isa Array{Float64, 3}
    @test size(inp_ss.Ietop, 1) == 2   # 2 beams
    @test size(inp_ss.Ietop, 2) == 1   # 1 time step
end

@testitem "find_input_file error cases" begin
    # No Ie_incoming file
    mktempdir() do emptydir
        @test_throws Exception find_input_file(emptydir)
    end

    # Multiple Ie_incoming files
    mktempdir() do twodir
        touch(joinpath(twodir, "Ie_incoming_a.mat"))
        touch(joinpath(twodir, "Ie_incoming_b.mat"))
        @test_throws Exception find_input_file(twodir)
    end
end
