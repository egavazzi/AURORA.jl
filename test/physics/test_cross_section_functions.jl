@testitem "Cross-section functions are order- and eltype-agnostic" begin
    E_grid = exp10.(range(log10(0.5), log10(1e5); length = 200))

    function_names = String[]
    for species in ("N2", "O2", "O")
        for state in AURORA.get_level_names(species)
            push!(function_names, "e_" * species * state)
        end
    end
    @test length(function_names) == 63

    # Deterministic, non-monotonic permutation: evens then odds.
    n = length(E_grid)
    p = vcat(2:2:n, 1:2:n)

    for fname in function_names
        f = getfield(AURORA, Symbol(fname))

        σ_sorted = f(E_grid)
        @test all(isfinite, σ_sorted)
        @test all(>=(0), σ_sorted)
        @test eltype(σ_sorted) <: AbstractFloat

        # Shuffled input matches sorted input permuted back.
        σ_shuffled = f(E_grid[p])
        @test σ_shuffled ≈ σ_sorted[p]

        # Int input matches Float64 input. Some functions are designed to throw for
        # energies below a fraction of an eV (see e_N2rot0_2), so keep rounded values >= 1.
        E_int = max.(1, round.(Int, E_grid))
        σ_int = f(E_int)
        σ_float = f(float.(E_int))
        @test eltype(σ_int) <: AbstractFloat
        @test σ_int == σ_float

        # A view over a reversed range preserves the input's axes.
        E_view = @view E_grid[end:-1:1]
        σ_view = f(E_view)
        @test axes(σ_view) == axes(E_view)
        @test σ_view ≈ reverse(σ_sorted)
    end
end

@testitem "Cross-section functions at splice and gate energies" begin
    # Energies where a function switches between two fits or interpolants, or where a gate
    # zeroes it. Each must be handled by exactly one branch.
    E_splice = [
        10.0^1.477121254719663,                 # e_N2rot0_2 PCHIP/linear splice
        6.867, 31.614,                          # e_O1D, e_O3p5P fit splices
        4.6469, 4.0781, 3.8671, 3.4811, 4.1979, # e_N2vib0_3…0_7 table/vib0_1-tail splices
        10.0, 25.0, 30.0, 200.0, 250.0, 300.0, 1000.0, 1e5,
        0.2888, 0.5742, 1.9, 1.9475, 1.68, 1.6801, 12.25, 12.255, # N2 gates
        16.1, 16.9,                                               # O2 gates
        4.17, 4.19, 10.73, 10.74, 13.6, 13.618,                   # O gates
    ]

    function_names = String[]
    for species in ("N2", "O2", "O")
        for state in AURORA.get_level_names(species)
            push!(function_names, "e_" * species * state)
        end
    end
    @test length(function_names) == 63

    for fname in function_names
        f = getfield(AURORA, Symbol(fname))
        for E in E_splice
            σ = f([E])
            @test length(σ) == 1
            @test isfinite(σ[1]) && σ[1] >= 0
            σ3 = f([E - 1e-6, E, E + 1e-6])
            @test σ[1] == σ3[2]
            if σ3[1] > 0 && σ3[3] > 0
                @test σ[1] > 0
            end
        end
    end
end
