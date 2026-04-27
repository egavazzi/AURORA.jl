@testmodule HelperFunc begin
    # Exact allowed-domain area for one primary bin.
    # For a given E_p, the allowed E_s range is [0, (E_p - E_threshold) / 2] and the
    # allowed E_d range is [(E_p - E_threshold) / 2, E_p - E_threshold].
    # Integrating their width (E_p - E_threshold) / 2 over the allowed E_p range [a, b]
    # gives
    #   \int_a^b (E_p - E_threshold) / 2 dE_p
    # = [(E_p - E_threshold)^2 / 4]_a^b
    # = ((b - E_threshold)^2 - (a - E_threshold)^2) / 4
    area_per_matrix(a, b, threshold) = ((b - threshold)^2 - (a - threshold)^2) / 4

    # Exact mean primary energy over the same allowed domain.
    # Higher E_p inside [a, b] contribute more because they allow a wider range of possible
    # energy splits, which is related to the excess energy (E_p - E_threshold).
    # As such, <E_p> is the weighted average using (E_p - E_threshold) as the weight.
    #    \int_a^b E_p * (E_p - E_threshold) dE_p / \int_a^b (E_p - E_threshold) dE_p
    function expected_primary_event_mean(a, b, threshold)
        numerator = ((b^3 - a^3) / 3) - threshold * ((b^2 - a^2) / 2)
        denominator = area_per_matrix(a, b, threshold) * 2
        return numerator / denominator
    end

    # Reflect a bin-integrated spectrum across E_excess / 2 using interval overlaps.
    # For non-uniform energy grid, each reflected bin mass is distributed to
    # the target bins in proportion to geometric overlap.
    function reflect_bin_weights(weights, E_edges, E_excess)
        reflected = zeros(eltype(weights), length(weights))

        for i_source in eachindex(weights)
            source_weight = weights[i_source]
            source_weight == 0 && continue

            source_left = E_edges[i_source]
            source_right = E_edges[i_source + 1]
            reflected_left = E_excess - source_right
            reflected_right = E_excess - source_left
            reflected_width = reflected_right - reflected_left
            reflected_width <= 0 && continue

            for i_target in eachindex(reflected)
                target_left = E_edges[i_target]
                target_right = E_edges[i_target + 1]
                overlap = min(reflected_right, target_right) - max(reflected_left, target_left)
                overlap <= 0 && continue

                reflected[i_target] += source_weight * overlap / reflected_width
            end
        end

        return reflected
    end


    # Refine a grid by splitting every bin at its midpoint, preserving nested bin edges.
    function midpoint_refined_edges(E_edges)
        midpoints = (E_edges[1:end-1] .+ E_edges[2:end]) ./ 2
        return sort!(vcat(E_edges, midpoints))
    end

    # Rebin bin-integrated weights onto a coarser target grid using interval overlaps.
    function rebin_weights(weights, source_edges, target_edges)
        rebinned = zeros(eltype(weights), length(target_edges) - 1)

        for i_source in eachindex(weights)
            source_weight = weights[i_source]
            source_weight == 0 && continue

            source_left = source_edges[i_source]
            source_right = source_edges[i_source + 1]
            source_width = source_right - source_left
            source_width <= 0 && continue

            for i_target in eachindex(rebinned)
                target_left = target_edges[i_target]
                target_right = target_edges[i_target + 1]
                overlap = min(source_right, target_right) - max(source_left, target_left)
                overlap <= 0 && continue

                rebinned[i_target] += source_weight * overlap / source_width
            end
        end

        return rebinned
    end

    # Fine-grid source bins that partition one coarse primary source bin.
    function child_bin_indices(coarse_edges, fine_edges, i_coarse)
        coarse_left = coarse_edges[i_coarse]
        coarse_right = coarse_edges[i_coarse + 1]
        indices = Int[]

        for i_fine in 1:(length(fine_edges) - 1)
            fine_left = fine_edges[i_fine]
            fine_right = fine_edges[i_fine + 1]
            if fine_left >= coarse_left && fine_right <= coarse_right
                push!(indices, i_fine)
            end
        end

        return indices
    end
end

@testmodule UniformSetup begin
    using AURORA

    THRESHOLD = 2.0
    E_EDGES = collect(0.0:0.125:8.0)
    E_CENTERS = (E_EDGES[1:end-1] .+ E_EDGES[2:end]) ./ 2
    # A constant-law where all energies have the same probability.
    SPEC = AURORA.CascadingSpec("TEST", [THRESHOLD], (E_s, E_p) -> 1.0)
    Q_PRIMARY, Q_SECONDARY, _, _ = AURORA.calculate_cascading_matrices(SPEC, E_EDGES; verbose = false)

    FIRST_ACTIVE_PRIMARY = findfirst(x -> x >= THRESHOLD, @view(E_EDGES[1:end-1]))
end

@testmodule NonUniformSetup begin
    using AURORA

    THRESHOLD = 16.0
    BASE_GRID = AURORA.EnergyGrid(1000.0)
    # Prepend 0 eV so we don't "lose" electrons that fall under the first bin edge of
    # the AURORA grid (typically around 2 eV). Important for conservation tests.
    lower_end = 0.0:0.1:BASE_GRID.E_edges[1]
    E_EDGES = vcat(lower_end[1:end-1], BASE_GRID.E_edges)
    E_CENTERS = (E_EDGES[1:end-1] .+ E_EDGES[2:end]) ./ 2
    # A constant-law where all energies have the same probability.
    SPEC = AURORA.CascadingSpec("TEST", [THRESHOLD], (E_s, E_p) -> 1.0)
    Q_PRIMARY, Q_SECONDARY, _, _ = AURORA.calculate_cascading_matrices(SPEC, E_EDGES; verbose = false)

    FIRST_ACTIVE_PRIMARY = findfirst(x -> x >= THRESHOLD, @view(E_EDGES[1:end-1]))
end




# Integration domains:
# Verifies, that the transfer matrices cover the full physically allowed domain.
@testitem "Cascading constant-law domain area" setup=[HelperFunc, UniformSetup] begin
    for i_primary in UniformSetup.FIRST_ACTIVE_PRIMARY:length(UniformSetup.E_CENTERS)
        a = UniformSetup.E_EDGES[i_primary]
        b = UniformSetup.E_EDGES[i_primary + 1]
        expected = HelperFunc.area_per_matrix(a, b, UniformSetup.THRESHOLD)

        # Constant law
        # -> unnormalized spectra are uniform over the allowed domain
        # -> their integrals are equal to the area of the domain.
        total_primary = sum(@view(UniformSetup.Q_PRIMARY[i_primary, :, 1]))
        total_secondary = sum(@view(UniformSetup.Q_SECONDARY[i_primary, :, 1]))

        @test isapprox(total_primary, expected; rtol = 1e-4)
        @test isapprox(total_secondary, expected; rtol = 1e-4)
    end
end
@testitem "Cascading nonuniform constant-law domain area" setup=[HelperFunc, NonUniformSetup] begin
    for i_primary in NonUniformSetup.FIRST_ACTIVE_PRIMARY:length(NonUniformSetup.E_CENTERS)
        a = NonUniformSetup.E_EDGES[i_primary]
        b = NonUniformSetup.E_EDGES[i_primary + 1]
        expected = HelperFunc.area_per_matrix(a, b, NonUniformSetup.THRESHOLD)

        total_primary = sum(@view(NonUniformSetup.Q_PRIMARY[i_primary, :, 1]))
        total_secondary = sum(@view(NonUniformSetup.Q_SECONDARY[i_primary, :, 1]))

        @test isapprox(total_primary, expected; rtol = 1e-4)
        @test isapprox(total_secondary, expected; rtol = 1e-4)
    end
end


# Particle conservation:
# Each ionization event should create one degraded and one secondary electron, so their
# unnormalized spectra should have equal integrals over the full domain.
@testitem "Cascading uniform particle conservation" setup=[UniformSetup] begin
    for i_primary in UniformSetup.FIRST_ACTIVE_PRIMARY:length(UniformSetup.E_CENTERS)
        total_primary = sum(@view(UniformSetup.Q_PRIMARY[i_primary, :, 1]))
        total_secondary = sum(@view(UniformSetup.Q_SECONDARY[i_primary, :, 1]))
        @test isapprox(total_primary, total_secondary; rtol = 1e-4, atol = 1e-12)
    end
end
@testitem "Cascading nonuniform particle conservation" setup=[NonUniformSetup] begin
    for i_primary in NonUniformSetup.FIRST_ACTIVE_PRIMARY:length(NonUniformSetup.E_CENTERS)
        total_primary = sum(@view(NonUniformSetup.Q_PRIMARY[i_primary, :, 1]))
        total_secondary = sum(@view(NonUniformSetup.Q_SECONDARY[i_primary, :, 1]))
        @test isapprox(total_primary, total_secondary; rtol = 1e-4, atol = 1e-12)
    end
end



# Energy conservation:
# The mean energy of the primary event should equal the mean degraded energy + mean
# secondary energy + the ionization threshold. The expected primary mean is calculated by
# weighting over the allowed (degraded+secondary) energy space.
@testitem "Cascading energy conservation" setup=[HelperFunc, UniformSetup] begin
    first_bin_width = UniformSetup.E_EDGES[2] - UniformSetup.E_EDGES[1]

    for i_primary in UniformSetup.FIRST_ACTIVE_PRIMARY:length(UniformSetup.E_CENTERS)
        a = UniformSetup.E_EDGES[i_primary]
        b = UniformSetup.E_EDGES[i_primary + 1]

        # Near threshold, the redistribution occupies too few bins for a bin-centre mean to
        # be a meaningful conservation diagnostic.
        (b - UniformSetup.THRESHOLD) <= 2 * first_bin_width && continue

        # wd / ws: unnormalized bin weights for degraded and secondary electrons
        # nd / ns: total weight
        wd = @view(UniformSetup.Q_PRIMARY[i_primary, :, 1])
        ws = @view(UniformSetup.Q_SECONDARY[i_primary, :, 1])
        nd = sum(wd)
        ns = sum(ws)
        nd == 0 && continue
        ns == 0 && continue

        E_d_mean = sum(wd .* UniformSetup.E_CENTERS) / nd
        E_s_mean = sum(ws .* UniformSetup.E_CENTERS) / ns
        E_p_event_mean = HelperFunc.expected_primary_event_mean(a, b, UniformSetup.THRESHOLD)

        @test isapprox(E_d_mean + E_s_mean + UniformSetup.THRESHOLD, E_p_event_mean;
                       rtol = 1e-3, atol = 0.01)
    end
end
@testitem "Cascading nonuniform energy conservation" setup=[HelperFunc, NonUniformSetup] begin
    first_bin_width = NonUniformSetup.E_EDGES[2]

    # for i_primary in (NonUniformSetup.FIRST_ACTIVE_PRIMARY + 5):length(NonUniformSetup.E_CENTERS)
    for i_primary in NonUniformSetup.FIRST_ACTIVE_PRIMARY:length(NonUniformSetup.E_CENTERS)
        a = NonUniformSetup.E_EDGES[i_primary]
        b = NonUniformSetup.E_EDGES[i_primary + 1]

        # Near threshold, the redistribution occupies too few bins for a bin-centre mean to
        # be a meaningful conservation diagnostic.
        (b - NonUniformSetup.THRESHOLD) <= 2 * first_bin_width && continue

        # wd / ws: unnormalized bin weights for degraded and secondary electrons
        # nd / ns: total weight
        wd = @view(NonUniformSetup.Q_PRIMARY[i_primary, :, 1])
        ws = @view(NonUniformSetup.Q_SECONDARY[i_primary, :, 1])
        nd = sum(wd)
        ns = sum(ws)
        nd == 0 && continue
        ns == 0 && continue

        E_d_mean = sum(wd .* NonUniformSetup.E_CENTERS) / nd
        E_s_mean = sum(ws .* NonUniformSetup.E_CENTERS) / ns
        E_p_event_mean = HelperFunc.expected_primary_event_mean(a, b, NonUniformSetup.THRESHOLD)

        @test isapprox(E_d_mean + E_s_mean + NonUniformSetup.THRESHOLD, E_p_event_mean;
                        rtol = 1e-3, atol = 0.01)
    end
end







# Symmetry of degraded and secondary spectra:
# For the constant-law case, the degraded and secondary spectra should be mirror images
# of each other around E_excess / 2. We test that by reflecting the secondary bin
# weights onto the degraded-energy axis using interval overlaps.
@testitem "Cascading degraded-secondary symmetry" setup=[HelperFunc, UniformSetup] begin
    # Check on one high-energy primary bin.
    i_primary = findfirst(i -> UniformSetup.E_EDGES[i] >= 6.0, eachindex(UniformSetup.E_CENTERS))

    # wd / ws: unnormalized bin weights for degraded and secondary electrons.
    # nd / ns: total weight (equal by particle conservation).
    # pd / ps: normalised probability distributions.
    wd = @view(UniformSetup.Q_PRIMARY[i_primary, :, 1])
    ws = @view(UniformSetup.Q_SECONDARY[i_primary, :, 1])
    nd = sum(wd)
    ns = sum(ws)

    @test nd > 0
    @test ns > 0

    pd = wd ./ nd
    ps = ws ./ ns
    E_excess = UniformSetup.E_CENTERS[i_primary] - UniformSetup.THRESHOLD
    reflected_ps = HelperFunc.reflect_bin_weights(ps, UniformSetup.E_EDGES, E_excess)

    # support = (pd .> 1e-12) .| (reflected_ps .> 1e-12)
    # @test sum(reflected_ps) > 0.95
    # @test isapprox(pd[support], reflected_ps[support]; atol = 0.01)

    @test sum(reflected_ps) > 0.95
    @test maximum(abs.(pd .- reflected_ps)) <= 0.01
end
@testitem "Cascading nonuniform degraded-secondary symmetry" setup=[HelperFunc, NonUniformSetup] begin
    # Check on one high-energy primary bin.
    i_primary = findfirst(i -> NonUniformSetup.E_EDGES[i] >= 600.0, eachindex(NonUniformSetup.E_CENTERS))

    # wd / ws: unnormalized bin weights for degraded and secondary electrons.
    # nd / ns: total weight (equal by particle conservation).
    # pd / ps: normalised probability distributions.
    wd = @view(NonUniformSetup.Q_PRIMARY[i_primary, :, 1])
    ws = @view(NonUniformSetup.Q_SECONDARY[i_primary, :, 1])
    nd = sum(wd)
    ns = sum(ws)

    @test nd > 0
    @test ns > 0

    pd = wd ./ nd
    ps = ws ./ ns
    E_excess = NonUniformSetup.E_CENTERS[i_primary] - NonUniformSetup.THRESHOLD
    reflected_ps = HelperFunc.reflect_bin_weights(ps, NonUniformSetup.E_EDGES, E_excess)

    # support = (pd .> 1e-12) .| (reflected_ps .> 1e-12)
    # @test sum(reflected_ps) > 0.95
    # @test isapprox(pd[support], reflected_ps[support]; atol = 0.03)

    @test sum(reflected_ps) > 0.95
    @test maximum(abs.(pd .- reflected_ps)) <= 0.01
end



# Integral boundaries:
# Verifies that bins outside the physical domain remain empty.
@testitem "Cascading domain edge behavior" setup=[UniformSetup] begin
    for i_primary in UniformSetup.FIRST_ACTIVE_PRIMARY:length(UniformSetup.E_CENTERS)
        a = UniformSetup.E_EDGES[i_primary]
        b = UniformSetup.E_EDGES[i_primary + 1]

        degraded_min = (a - UniformSetup.THRESHOLD) / 2
        degraded_max = b - UniformSetup.THRESHOLD
        secondary_max = (b - UniformSetup.THRESHOLD) / 2

        for i_bin in eachindex(UniformSetup.E_CENTERS)
            E_bin_min = UniformSetup.E_EDGES[i_bin]
            E_bin_max = UniformSetup.E_EDGES[i_bin + 1]

            primary_outside_domain = (E_bin_max <= degraded_min) || (E_bin_min >= degraded_max)
            secondary_outside_domain = E_bin_min >= secondary_max

            if primary_outside_domain
                @test UniformSetup.Q_PRIMARY[i_primary, i_bin, 1] ≈ 0.0
            end
            if secondary_outside_domain
                @test UniformSetup.Q_SECONDARY[i_primary, i_bin, 1] ≈ 0.0
            end
        end
    end
end
@testitem "Cascading nonuniform domain edge behavior" setup=[NonUniformSetup] begin
    for i_primary in NonUniformSetup.FIRST_ACTIVE_PRIMARY:length(NonUniformSetup.E_CENTERS)
        a = NonUniformSetup.E_EDGES[i_primary]
        b = NonUniformSetup.E_EDGES[i_primary + 1]

        degraded_min = (a - NonUniformSetup.THRESHOLD) / 2
        degraded_max = b - NonUniformSetup.THRESHOLD
        secondary_max = (b - NonUniformSetup.THRESHOLD) / 2

        for i_bin in eachindex(NonUniformSetup.E_CENTERS)
            E_bin_min = NonUniformSetup.E_EDGES[i_bin]
            E_bin_max = NonUniformSetup.E_EDGES[i_bin + 1]

            primary_outside_domain = (E_bin_max <= degraded_min) || (E_bin_min >= degraded_max)
            secondary_outside_domain = E_bin_min >= secondary_max

            if primary_outside_domain
                @test NonUniformSetup.Q_PRIMARY[i_primary, i_bin, 1] ≈ 0.0
            end
            if secondary_outside_domain
                @test NonUniformSetup.Q_SECONDARY[i_primary, i_bin, 1] ≈ 0.0
            end
        end
    end
end




# Grid convergence:
# Verifies that refining the grid gives consistent bin-integrated spectra. For each coarse
# primary bin, we sum the fine-grid spectra of its child bins and compare them to the
# coarse-grid spectrum.
@testitem "Cascading uniform grid convergence" setup=[HelperFunc, UniformSetup] begin
    using AURORA

    coarse_edges = UniformSetup.E_EDGES
    fine_edges = HelperFunc.midpoint_refined_edges(coarse_edges)

    Qp_coarse, Qs_coarse, _, _ = AURORA.calculate_cascading_matrices(UniformSetup.SPEC, coarse_edges; verbose = false)
    Qp_fine, Qs_fine, _, _ = AURORA.calculate_cascading_matrices(UniformSetup.SPEC, fine_edges; verbose = false)

    first_active = findfirst(x -> x >= UniformSetup.THRESHOLD, @view(coarse_edges[1:end-1]))

    for i_primary in first_active:(length(coarse_edges) - 1)
        child_indices = HelperFunc.child_bin_indices(coarse_edges, fine_edges, i_primary)
        isempty(child_indices) && continue

        coarse_primary = @view(Qp_coarse[i_primary, :, 1])
        coarse_secondary = @view(Qs_coarse[i_primary, :, 1])

        fine_primary = sum(Qp_fine[child_indices, :, 1]; dims = 1)[:]
        fine_secondary = sum(Qs_fine[child_indices, :, 1]; dims = 1)[:]

        rebinned_primary = HelperFunc.rebin_weights(fine_primary, fine_edges, coarse_edges)
        rebinned_secondary = HelperFunc.rebin_weights(fine_secondary, fine_edges, coarse_edges)

        @test isapprox(coarse_primary, rebinned_primary; rtol = 1e-4)
        @test isapprox(coarse_secondary, rebinned_secondary; rtol = 1e-4)
    end
end
@testitem "Cascading nonuniform grid convergence" setup=[HelperFunc, NonUniformSetup] begin
    using AURORA

    coarse_edges = NonUniformSetup.E_EDGES
    fine_edges = HelperFunc.midpoint_refined_edges(coarse_edges)

    Qp_coarse, Qs_coarse, _, _ = AURORA.calculate_cascading_matrices(NonUniformSetup.SPEC, coarse_edges; verbose = false)
    Qp_fine, Qs_fine, _, _ = AURORA.calculate_cascading_matrices(NonUniformSetup.SPEC, fine_edges; verbose = false)

    first_active = findfirst(x -> x >= NonUniformSetup.THRESHOLD, @view(coarse_edges[1:end-1]))

    for i_primary in first_active:(length(coarse_edges) - 1)
        child_indices = HelperFunc.child_bin_indices(coarse_edges, fine_edges, i_primary)
        isempty(child_indices) && continue

        coarse_primary = @view(Qp_coarse[i_primary, :, 1])
        coarse_secondary = @view(Qs_coarse[i_primary, :, 1])

        fine_primary = sum(Qp_fine[child_indices, :, 1]; dims = 1)[:]
        fine_secondary = sum(Qs_fine[child_indices, :, 1]; dims = 1)[:]

        rebinned_primary = HelperFunc.rebin_weights(fine_primary, fine_edges, coarse_edges)
        rebinned_secondary = HelperFunc.rebin_weights(fine_secondary, fine_edges, coarse_edges)

        @test isapprox(coarse_primary, rebinned_primary; rtol = 1e-4)
        @test isapprox(coarse_secondary, rebinned_secondary; rtol = 1e-4)
    end
end
