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

    raw_t = Float64[]
    let previous_t_last = nothing
        for file in AURORA.list_result_files(SharedSimResults.td_dir)
            _, t_local = AURORA.read_results(file; previous_t_last=previous_t_last)
            append!(raw_t, collect(t_local))
            previous_t_last = isempty(t_local) ? previous_t_last : t_local[end]
        end
    end
    @test vol_td.t == raw_t
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

    raw_t = Float64[]
    let previous_t_last = nothing
        for file in AURORA.list_result_files(SharedSimResults.td_dir)
            _, t_local = AURORA.read_results(file; previous_t_last=previous_t_last)
            append!(raw_t, collect(t_local))
            previous_t_last = isempty(t_local) ? previous_t_last : t_local[end]
        end
    end
    @test col_td.t == raw_t
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

@testitem "read_results drops legacy duplicate boundaries" begin
    using MAT

    mktempdir() do savedir
        file1 = joinpath(savedir, "IeFlickering-01.mat")
        io = matopen(file1, "w")
            write(io, "Ie_ztE", reshape(collect(1.0:12.0), 2, 3, 2))
            write(io, "t_run", [0.0, 0.01, 0.02])
        close(io)

        file2 = joinpath(savedir, "IeFlickering-02.mat")
        io = matopen(file2, "w")
            write(io, "Ie_ztE", reshape(collect(13.0:24.0), 2, 3, 2))
            write(io, "t_run", [0.02, 0.03, 0.04])
        close(io)

        files = AURORA.list_result_files(savedir)
        Ie1, t1 = AURORA.read_results(files[1])
        Ie2, t2 = AURORA.read_results(files[2]; previous_t_last=t1[end])

        @test collect(t1) == [0.0, 0.01, 0.02]
        @test collect(t2) == [0.03, 0.04]
        @test size(Ie1, 2) == 3
        @test size(Ie2, 2) == 2
    end
end
