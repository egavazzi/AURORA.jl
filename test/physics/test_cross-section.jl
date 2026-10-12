@testitem "Built-in channel tables build the species matrices" begin
    energy_grid = AURORA.EnergyGrid(500)

    for species in (:N2, :O2, :O)
        d = AURORA.default_collisions(species)
        channels = d.channels
        σ = AURORA.channel_cross_sections(d.elastic_cross_section, channels,
                                          energy_grid.E_centers)
        levels = AURORA.channel_excitation_levels(channels)

        # Row 1 is elastic, row i+1 is channels[i]
        @test size(σ) == (length(channels) + 1, energy_grid.n)
        @test size(levels) == (length(channels) + 1, 2)
        @test σ[1, :] == d.elastic_cross_section(energy_grid.E_centers)
        @test levels[1, :] == [0.0, 0.0]
        @test levels[2:end, 1] == [c.energy_loss for c in channels]
        @test levels[2:end, 2] == [c.n_secondaries for c in channels]
        @test all(c -> c.source == "built-in", channels)
    end
    @test_throws "No built-in collision channels for Ar" AURORA.default_collisions(:Ar)

    # Hard-coded values: any edit of the tables shows up here
    N2 = AURORA.default_collisions(:N2).channels
    @test AURORA.channel_names(N2)[1:4] == ["rot0_2", "rot0_4", "rot0_6", "rot0_8"]
    @test AURORA.channel_names(N2)[end-4:end] ==
          ["ionx2sgp", "iona2pu", "ionb2sup", "dion", "ddion"]
    @test [c.energy_loss for c in N2] == [
        0.00148010556, 0.004933884, 0.01036181244, 0.01776464064,
        0.2888, 0.5742, 0.8559, 1.1342, 1.4088, 1.6801, 1.9475,
        6.1688, 7.3532, 7.3622, 8.1647, 8.3987, 8.8895, 11.875, 12.255,
        8.5489, 11.032, 12.85, 12.94, 12.08, 12.85, 12.75, 12.8, 13.15, 13.1,
        20.6, 15.581, 16.73, 18.75, 24.0, 42.0]
    @test [c.n_secondaries for c in N2] == [zeros(Int, 30); 1; 1; 1; 1; 2]

    O2 = AURORA.default_collisions(:O2).channels
    O = AURORA.default_collisions(:O).channels
    @test [c.energy_loss for c in AURORA.ionizing_channels(O2)] ==
          [12.072, 16.1, 16.9, 18.2, 18.9, 32.51]
    @test [c.energy_loss for c in AURORA.ionizing_channels(O)] == [13.618, 16.9, 18.6, 28.5]

    # A few cross-section values (m²) at a fixed energy
    E = [100.0]
    @test AURORA.channel(N2, "ionx2sgp").cross_section(E)[1] ≈ 8.68168140019684e-21 rtol=1e-12
    @test AURORA.default_collisions(:N2).elastic_cross_section(E)[1] ≈
          5.513124309180136e-20 rtol=1e-12
end

@testitem "Cascading specs are derived from the channel tables" begin
    spec = AURORA.default_cascading_spec(:N2)
    @test spec.name == "N2"
    @test spec.ionization_thresholds == [15.581, 16.73, 18.75, 24.0, 42.0]
    @test spec.n_secondaries == [1, 1, 1, 1, 2]

    # Channels sharing an energy loss must agree on the number of secondaries
    σ = AURORA.default_collisions(:N2).elastic_cross_section
    law = AURORA.default_collisions(:N2).secondary_law
    clashing = [AURORA.CollisionChannel("a", σ, 20.0, 1),
                AURORA.CollisionChannel("b", σ, 20.0, 2)]
    @test_throws "must eject the same number of secondaries" AURORA.cascading_spec_from_channels(
        "Clash", law, clashing)

    # Sharing an energy loss is fine when the secondary counts agree: one threshold
    agreeing = [AURORA.CollisionChannel("a", σ, 20.0, 1),
                AURORA.CollisionChannel("b", σ, 20.0, 1)]
    @test AURORA.cascading_spec_from_channels("Agree", law,
                                              agreeing).ionization_thresholds == [20.0]

    # Without ionizing channels a spec needs no law; with one it does
    excitation = [AURORA.CollisionChannel("a", σ, 20.0, 0)]
    spec = AURORA.cascading_spec_from_channels("Exc", nothing, excitation)
    @test isempty(spec.ionization_thresholds)
    @test_throws "has ionizing channels (b) but no secondary_law" (
        AURORA.cascading_spec_from_channels("Ion", nothing, [excitation; agreeing[2]]))
end

@testitem "CollisionChannel validates its inputs" begin
    σ = AURORA.default_collisions(:N2).elastic_cross_section

    @test_throws "bare anonymous function" AURORA.CollisionChannel("x", E -> E, 1.0, 0)
    @test_throws "must be finite and non-negative" AURORA.CollisionChannel("x", σ, -1.0, 0)
    @test_throws "0 (excitation)" AURORA.CollisionChannel("x", σ, 1.0, 3)
    @test_throws "positive energy loss" AURORA.CollisionChannel("x", σ, 0.0, 1)

    c = AURORA.CollisionChannel("x", σ, 1.0, 0; source = "made up")
    @test c.source == "made up"
    c2 = AURORA.CollisionChannel(c; energy_loss = 14.0, n_secondaries = 1)
    @test c2.name == "x" && c2.energy_loss == 14.0 && c2.n_secondaries == 1
    @test c2.cross_section === c.cross_section && c2.source == "made up"

    @test_throws KeyError AURORA.channel([c], "absent")
    @test AURORA.channel([c], "x") === c
    @test_throws "Multiple channels are named x" AURORA.channel([c, c2], "x")
    @test AURORA.channel_names([c, c2]) == ["x", "x"]
    @test AURORA.ionizing_channels([c, c2]) == [c2]
end

@testitem "replace_channel! swaps one channel in place" begin
    channels = AURORA.default_collisions(:N2).channels
    old = AURORA.channel(channels, "a3sup")
    new = AURORA.replace_channel!(channels, "a3sup"; energy_loss = 6.2, source = "edited")
    @test AURORA.channel(channels, "a3sup") === new
    @test new.energy_loss == 6.2 && new.source == "edited"
    @test new.cross_section === old.cross_section && new.n_secondaries == 0
    @test AURORA.channel_names(channels) ==
          AURORA.channel_names(AURORA.default_collisions(:N2).channels)
    @test_throws KeyError AURORA.replace_channel!(channels, "absent"; energy_loss = 1.0)

    msis_file = find_msis_file(; verbose = false)
    sp = AURORA.N2Species(read_msis_file(msis_file)[:N2])
    new = AURORA.replace_channel!(sp, "ddion"; n_secondaries = 1)
    @test AURORA.channel(sp, "ddion") === new
    @test new.n_secondaries == 1
    @test_throws KeyError AURORA.replace_channel!(sp, "absent"; energy_loss = 1.0)
end

@testitem "Species fields reject laws that cannot be saved" begin
    msis_file = find_msis_file(; verbose = false)
    sp = AURORA.N2Species(read_msis_file(msis_file)[:N2])

    @test_throws "bare anonymous function" sp.secondary_law = (E_s, E_p) -> E_s
    @test_throws "bare anonymous function" sp.elastic_cross_section = E -> E
    @test_throws "bare anonymous function" sp.density_source = h -> h
    @test_throws "bare anonymous function" sp.phase_fcn_generator = (θ, E) -> θ

    # The rejected assignments left the species untouched
    @test sp.elastic_cross_section === AURORA.default_collisions(:N2).elastic_cross_section

    # Named functions, functors and @law are accepted
    sp.elastic_cross_section = AURORA.e_O2elastic
    @test sp.elastic_cross_section === AURORA.e_O2elastic
    sp.secondary_law = @law (E_s, E_p) -> 1.0 / (10.0^2 + E_s^2)
    @test sp.secondary_law isa ExprLaw
end
