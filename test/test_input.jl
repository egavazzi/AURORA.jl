using AURORA

@testset "Maxwellian with LET" begin
    ## Setting parameters
    altitude_lims = [100, 600];     # (km) altitude limits of the ionosphere
    θ_lims = 180:-10:0              # (°) angle-limits for the electron beams
    E_max = 100000;                 # (eV) upper limit to the energy grid
    msis_file = find_nrlmsis_file();
    iri_file = find_iri_file();
    h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
    μ_scatterings = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file);

    # Check that IeE_top is respected for a Maxwellian without LET (with a tolerance of 0.1%)
    Ie_top = Ie_with_LET(100, 1e5, E, dE, μ_center, μ_scatterings.BeamWeight, 1, low_energy_tail = false)
    IeE_top = sum(Ie_top .* reshape(E .+ dE / 2, (1, 1, :)))
    @test isapprox(IeE_top, 1e5, rtol = 0.001)

    # Check that varying the beams does not change IeE_top
    Ie_top_LET_1 = Ie_with_LET(100, 1e5, E, dE, μ_center, μ_scatterings.BeamWeight, 1, low_energy_tail = true)
    Ie_top_LET_2 = Ie_with_LET(100, 1e5, E, dE, μ_center, μ_scatterings.BeamWeight, 1:2, low_energy_tail = true)
    Ie_top_LET_3 = Ie_with_LET(100, 1e5, E, dE, μ_center, μ_scatterings.BeamWeight, [2, 5], low_energy_tail = true)
    IeE_top_1 = sum(Ie_top_LET_1 .* reshape(E .+ dE / 2, (1, 1, :)))
    IeE_top_2 = sum(Ie_top_LET_2 .* reshape(E .+ dE / 2, (1, 1, :)))
    IeE_top_3 = sum(Ie_top_LET_3 .* reshape(E .+ dE / 2, (1, 1, :)))

    @test IeE_top_1 ≈ IeE_top_2 ≈ IeE_top_3
end








## Code to reproduce Fig 4. of Meier et al. 1989
# using GLMakie
# # Setting parameters
# altitude_lims = [100, 600];     # (km) altitude limits of the ionosphere
# θ_lims = 180:-10:0              # (°) angle-limits for the electron beams
# E_max = 100000;                 # (eV) upper limit to the energy grid
# msis_file = find_nrlmsis_file();
# iri_file = find_iri_file();
# h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
# μ_scatterings = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file);
# # Plotting
# f = Figure()
# ax = Axis(f[1, 1], xscale = log10, yscale = log10)
# for (i, E₀) in enumerate([100, 250, 500, 1e3, 2e3, 4e3, 8e3])
#     Ie_top = Ie_with_LET(E₀, 6.242e15, E, dE, μ_center, μ_scatterings.BeamWeight, 1, low_energy_tail = false)
#     Ie_top_withtail = Ie_with_LET(E₀, 6.242e15, E, dE, μ_center, μ_scatterings.BeamWeight, 1, low_energy_tail = true)
#     lines!(E, Ie_top[1, 1, :] ./ dE, color = Cycled(i))
#     lines!(E, Ie_top_withtail[1, 1, :] ./ dE, linestyle = :dash, color = Cycled(i))
# end
# ylims!(10 * 1e4, 1e8 * 1e4)
# xlims!(10, nothing)
# display(f)
