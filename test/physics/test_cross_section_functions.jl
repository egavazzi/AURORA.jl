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
        @test σ_shuffled == σ_sorted[p]

        # Int input matches Float64 input. `e_N2rot0_2` throws below ~0.03 eV, so keep
        # rounded values >= 1.
        E_int = max.(1, round.(Int, E_grid))
        σ_int = f(E_int)
        σ_float = f(float.(E_int))
        @test eltype(σ_int) <: AbstractFloat
        @test σ_int == σ_float

        # A view over a reversed range preserves the input's axes.
        E_view = @view E_grid[end:-1:1]
        σ_view = f(E_view)
        @test axes(σ_view) == axes(E_view)
        @test σ_view == reverse(σ_sorted)
    end
end

@testitem "Cross sections are zero below their energy loss" begin
    E_centers = [1.0, 2.0, 3.0, 4.0]
    σ = [1.0, 2.0, 0.0, 4.0]
    @test_logs (:warn, r"non-zero cross section below its energy loss of 3.5 eV") AURORA.zero_below_energy_loss!(
        σ, E_centers, 3.5, "test channel")
    @test σ == [0.0, 0.0, 0.0, 4.0]

    # Nothing to zero: no warning, values unchanged.
    σ = [0.0, 0.0, 1.0, 2.0]
    @test_logs AURORA.zero_below_energy_loss!(σ, E_centers, 2.5, "another test channel")
    @test σ == [0.0, 0.0, 1.0, 2.0]
end
