@testitem "Cross-section functions are order- and eltype-agnostic" begin
    using OffsetArrays
    E_grid = exp10.(range(log10(0.5), log10(1e5); length = 200))

    cross_section_functions = Any[]
    for species in (:N2, :O2, :O)
        push!(cross_section_functions, AURORA.default_elastic_cross_section(species))
        for c in AURORA.default_channels(species)
            push!(cross_section_functions, c.cross_section)
        end
    end
    @test length(cross_section_functions) == 63

    # Deterministic, non-monotonic permutation: evens then odds.
    n = length(E_grid)
    p = vcat(2:2:n, 1:2:n)

    for f in cross_section_functions
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

        # An offset vector gets an offset result with the same axes.
        E_offset = OffsetArray(E_grid[p], -5)
        σ_offset = f(E_offset)
        @test axes(σ_offset) == axes(E_offset)
        @test OffsetArrays.no_offset_view(σ_offset) == σ_sorted[p]
    end
end
