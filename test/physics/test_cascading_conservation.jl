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
    SPEC = AURORA.CascadingSpec("TEST", [THRESHOLD], @law((E_s, E_p) -> 1.0))
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
    SPEC = AURORA.CascadingSpec("TEST", [THRESHOLD], @law((E_s, E_p) -> 1.0))
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

        # Coarse and rebinned-fine entries each carry quadrature error at rtol = 1e-4
        # (see fill_single_ionization_bin!), so they agree only to a few 1e-4.
        @test isapprox(coarse_primary, rebinned_primary; rtol = 1e-3)
        @test isapprox(coarse_secondary, rebinned_secondary; rtol = 1e-3)
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

        # Coarse and rebinned-fine entries each carry quadrature error at rtol = 1e-4
        # (see fill_single_ionization_bin!), so they agree only to a few 1e-4
        @test isapprox(coarse_primary, rebinned_primary; rtol = 1e-3)
        @test isapprox(coarse_secondary, rebinned_secondary; rtol = 1e-3)
    end
end


# Real-law energy conservation (to catch cases that the constant-law tests might miss)
@testitem "Cascading real-law energy conservation (floored grid)" begin
    using AURORA

    spec = AURORA.DefaultCascadingSpecN2()             # real law 1 / (11.4² + E_s²)
    I = spec.ionization_thresholds[1]                  # 15.581 eV
    eg = AURORA.EnergyGrid(3000.0)                     # standard grid, floor ≈ 2 eV
    Ec = eg.E_centers
    nE = length(Ec)

    # Force a fresh recompute on this grid
    cache = AURORA.SpeciesCascadingCache(spec)
    policy = AURORA.CachePolicy(force_recompute = true, save_cache = false)
    AURORA.load_or_compute_cascading!(cache, eg; verbose = false, policy)

    # A single ionizing channel (threshold I, one secondary) with unit cross-section. Level 1 is
    # the (skipped) elastic placeholder, matching the real E_levels layout.
    E_levels = [0.0 0.0; I 1.0]
    σ = zeros(2, nE); σ[2, :] .= 1.0

    # Per-event placed-energy ratio for every primary bin above threshold. (push! mutates
    # `ratios` in place — avoids the soft-scope pitfall of `+=` on an outer var inside a loop.)
    ratios = Float64[]
    for iE in eachindex(Ec)
        Ec[iE] - I < 10 && continue                    # skip near-threshold (too few bins to be meaningful)
        sec  = zeros(nE)
        prim = zeros(nE)
        AURORA.compute_ionization_spectra!(sec, prim, σ, E_levels, cache, iE)
        sum(prim) <= 0 && continue
        # Energy carried away by the outgoing electrons (degraded primary + secondaries) per event.
        placed = sum(Ec .* (sec .+ prim)) / σ[2, iE]
        push!(ratios, placed / (Ec[iE] - I))
    end

    @test length(ratios) > 50                          # ensure we actually exercised the spectrum
    # The cascade must never create energy: outgoing electrons carry ≤ the available excess.
    @test maximum(ratios) <= 1.005
end


# Double-ionization energy conservation
@testitem "Cascading double-ionization energy conservation (floored grid)" begin
    using AURORA

    spec = AURORA.DefaultCascadingSpecN2()             # real law 1 / (11.4² + E_s²)
    @test spec.n_secondaries[end] == 2                 # sanity: last channel is double ionization
    I = spec.ionization_thresholds[end]                # 42.0 eV double-ionization threshold
    eg = AURORA.EnergyGrid(3000.0)                     # standard grid, floor ≈ 2 eV
    Ec = eg.E_centers
    nE = length(Ec)

    # Force a fresh recompute on this grid
    cache = AURORA.SpeciesCascadingCache(spec)
    policy = AURORA.CachePolicy(force_recompute = true, save_cache = false)
    AURORA.load_or_compute_cascading!(cache, eg; verbose = false, policy)

    # A single DOUBLE-ionizing channel: threshold I ejects two secondaries (column 2 == 2). Level 1
    # is the (skipped) elastic placeholder, matching the real E_levels layout.
    E_levels = [0.0 0.0; I 2.0]
    σ = zeros(2, nE); σ[2, :] .= 1.0

    ratios = Float64[]
    for iE in eachindex(Ec)
        Ec[iE] - I < 10 && continue                    # skip near-threshold (too few bins to be meaningful)
        sec  = zeros(nE)
        prim = zeros(nE)
        AURORA.compute_ionization_spectra!(sec, prim, σ, E_levels, cache, iE)
        sum(prim) <= 0 && continue
        # Energy carried away by the three outgoing electrons (degraded primary + 2 secondaries).
        placed = sum(Ec .* (sec .+ prim)) / σ[2, iE]
        push!(ratios, placed / (Ec[iE] - I))
    end

    @test length(ratios) > 50
    # Never create energy: degraded primary + two secondaries carry ≤ the available excess.
    @test maximum(ratios) <= 1.005
    # And well above threshold (secondaries mostly on-grid) the cascade places ~all of the excess,
    # confirming the two secondaries are actually deposited rather than dropped.
    @test maximum(ratios) >= 0.95
end

# A law rebuilt from its source (as happens when a saved model is reloaded) is evaluated in a
# newer world age than the frame that goes on to use it. Building the matrices from inside a
# single function reproduces that situation; at top level each statement starts a new world,
# so the failure only shows up here. The result must match a plain functor law exactly.
@testitem "Cascading matrices from a law rebuilt in a newer world age" setup=[UniformSetup] begin
    using AURORA
    struct FlatLaw end
    (::FlatLaw)(E_s, E_p) = 1.0

    function build_from_source()
        law = AURORA.ExprLaw("(E_s, E_p) -> 1.0")
        spec = AURORA.CascadingSpec("TEST", [UniformSetup.THRESHOLD], law)
        return AURORA.calculate_cascading_matrices(spec, UniformSetup.E_EDGES; verbose = false)
    end
    Qp_src, Qs_src, _, _ = build_from_source()

    spec_functor = AURORA.CascadingSpec("TEST", [UniformSetup.THRESHOLD], FlatLaw())
    Qp_fn, Qs_fn, _, _ = AURORA.calculate_cascading_matrices(spec_functor, UniformSetup.E_EDGES; verbose = false)

    @test Qp_src == Qp_fn
    @test Qs_src == Qs_fn
    @test Qp_src == UniformSetup.Q_PRIMARY
    @test Qs_src == UniformSetup.Q_SECONDARY

# The production double-ionization path (numerical-CDF + fixed Gauss-Legendre rules,
# `fill_double_ionization_bin_cdf!`) has no built-in error estimate, so validate it here
# against the adaptive-cubature reference (`fill_double_ionization_bin!`, rtol = 1e-3),
# which has a controlled error for any secondary law. Two laws of different regularity are
# checked: the smooth Lorentzian-like N2/O2 shape and the atomic-O law, whose
# piecewise-linear E_p-dependence has kinks at its parameter knots. Tolerances cover both
# methods' quadrature errors: on rows carrying real weight the row sums agree to a few
# 1e-3 and individual entries to ~1e-2; rows just above the threshold carry ~1e-5 of the
# weight and are allowed a looser bound (the fixed primary-bin rule does not resolve the
# kinematic-clamp kinks there).
@testitem "Cascading double-ionization CDF matches adaptive reference" begin
    using AURORA
    using HCubature: hcubature_buffer

    smooth_law = AURORA.@law (E_s, E_p) -> 1.0 / (11.4^2 + E_s^2)
    o_law = AURORA.OSecondaryLaw([100.0, 200, 500, 1000, 2000],
                                 [12.6, 13.7, 14.1, 14.0, 13.7],
                                 [7.18, 4.97, 2.75, 1.69, 1.02] .* 1e-22)

    E_edges, _, _ = AURORA.make_energy_grid(1500.0)
    E_left = E_edges[1:end-1]
    n_E = length(E_left)
    # Grid for the event-count consistency check: subdivide [0, grid floor) into thirds so
    # no secondary energy falls below the lowest edge.
    zero_floor_edges = vcat(range(0.0, E_edges[1]; length = 4)[1:end-1], E_edges)

    for (law, threshold) in ((smooth_law, 42.0), (o_law, 28.5))
        spec = AURORA.CascadingSpec("TEST", [threshold], law; n_secondaries = [2])
        i_first = searchsortedfirst(E_left, threshold)

        # Production path
        Qp, Qs, _, _ = AURORA.calculate_cascading_matrices(spec, E_edges; verbose = false)

        # Adaptive reference, on every active row
        P = zeros(n_E, n_E, 1)
        S = zeros(n_E, n_E, 1)
        pbuf = hcubature_buffer(AURORA.DoublePrimaryCascadingIntegrand(0.0, 1.0, 0.0, law),
                                (0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
        sbuf = hcubature_buffer(AURORA.DoubleSecondaryCascadingIntegrand(0.0, 1.0, 0.0, law),
                                (0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
        for i_p in i_first:n_E
            AURORA.fill_double_ionization_bin!(P, S, E_edges, E_left, threshold, i_p, 1,
                                               law, pbuf, sbuf)
        end

        for (ref, cdf) in ((P, Qp), (S, Qs))
            row_ref = vec(sum(ref[:, :, 1]; dims = 2))
            row_cdf = vec(sum(cdf[:, :, 1]; dims = 2))
            w_max = maximum(row_ref)
            for i in i_first:n_E
                row_ref[i] > 0 || continue
                rel = abs(row_cdf[i] - row_ref[i]) / max(row_ref[i], row_cdf[i])
                if row_ref[i] > 0.01 * w_max
                    @test rel < 2e-2   # rows that matter: both methods within a few 1e-3
                    # The solver consumes the spectra bin by bin, so the row's shape matters
                    # too: check entries carrying at least 1 % of the row's largest entry.
                    # The bins clipped by the kinematic boundaries and the rows just above
                    # the threshold see the clamp kinks the fixed rules do not resolve
                    # (measured up to ~1.2e-1); everywhere else both methods agree to ~1e-2.
                    e_max = maximum(@view ref[i, :, 1])
                    jlo, jhi = extrema(findall(>(0), @view ref[i, :, 1]))
                    for j in 1:n_E
                        ref[i, j, 1] > 0.01 * e_max || continue
                        rel_entry = abs(cdf[i, j, 1] - ref[i, j, 1]) /
                                    max(ref[i, j, 1], cdf[i, j, 1])
                        if j <= jlo + 1 || j >= jhi - 1 || E_left[i] < 4 * threshold
                            @test rel_entry < 2e-1
                        else
                            @test rel_entry < 3e-2
                        end
                    end
                else
                    @test rel < 2e-1   # near-threshold rows: negligible weight, loose bound
                end
            end
        end

        # Internal consistency: both marginals count the same double-ionization events
        # (Z₂), so their row sums must agree wherever all three outgoing electrons land
        # on-grid. On the standard grid the secondary marginal loses the below-floor part
        # of the law by design, so run this check on the zero-floor grid.
        Qp0, Qs0, _, _ = AURORA.calculate_cascading_matrices(spec, zero_floor_edges;
                                                             verbose = false)
        E_left0 = zero_floor_edges[1:end-1]
        row_p = vec(sum(Qp0[:, :, 1]; dims = 2))
        row_s = vec(sum(Qs0[:, :, 1]; dims = 2))
        i_hi = searchsortedfirst(E_left0, 10 * threshold)
        for i in i_hi:length(E_left0)
            row_p[i] > 0 || continue
            @test isapprox(row_p[i], row_s[i]; rtol = 5e-3)
        end
    end
end
