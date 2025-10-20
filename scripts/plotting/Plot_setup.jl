using AURORA
using CairoMakie
CairoMakie.activate!()

altitude_lims = [100, 600];     # (km) altitude limits of the ionosphere
θ_lims = 180:-10:0          # (°) angle-limits for the electron beams
E_max = 3000;               # (eV) upper limit to the energy grid

msis_file = find_msis_file(year = 2005, month = 10, day = 8, hour = 22, minute = 0,
                              lat = 70, lon = 19, height = 85:1:700);
iri_file = find_iri_file(year = 2005, month = 10, day = 8, hour = 22, minute = 0, lat = 70,
                         lon = 19, height = 85:1:700);

h_atm, ne, Te, Tn, E, dE, n_neutrals, E_levels_neutrals, σ_neutrals, μ_lims, μ_center,
μ_scatterings = setup(altitude_lims, θ_lims, E_max, msis_file, iri_file);

## Plot densities
f1 = Figure()
ax11 = Axis(f1[1,1], xlabel = "(m⁻³)", ylabel = "height (km)",
          title = "MSIS neutral density", xscale = log10, xticks = LogTicks(LinearTicks(9)))
lines!(n_neutrals.nN2, h_atm/1e3, linewidth = 2, label = "N₂")
lines!(n_neutrals.nO2, h_atm/1e3, linestyle = :dash, linewidth = 2, color = :red, label = "O₂")
lines!(n_neutrals.nO, h_atm/1e3, linestyle = :dashdot, linewidth = 2, label = "O")
xlims!(1e12, 1e20)
ylims!(85, 700)
axislegend("Densities", position = :rt)
ax12 = Axis(f1[1,2], xlabel = "(m⁻³)", ylabel = "height (km)",
        title = "electron density", xscale = log10)
lines!(ne, h_atm/1e3, linewidth = 2)
xlims!(1e8, 1e13)
ylims!(85, 700)
display(f1)

## Plot excitation thresholds
f2 = Figure()
ax21 = Axis(f2[1,1], ylabel = "Excitation energy (eV)", title = "N₂")
ax22 = Axis(f2[1,2], title = "O₂")
ax23 = Axis(f2[1,3], title = "O")
for i in eachindex(E_levels_neutrals.N2_levels[:,1])
lines!(ax21, [0, 1], E_levels_neutrals.N2_levels[i, 1].*[1, 1], linewidth = 2)
end
for i in eachindex(E_levels_neutrals.O2_levels[:,1])
lines!(ax22, [0, 1], E_levels_neutrals.O2_levels[i, 1].*[1, 1], linewidth = 2)
end
for i in eachindex(E_levels_neutrals.O_levels[:,1])
lines!(ax23, [0, 1], E_levels_neutrals.O_levels[i, 1].*[1, 1], linewidth = 2)
end
[xlims!(X, 0, 1) for X in [ax21, ax22, ax23]]
display(f2)

## Plot the energy grid
f3 = Figure()
ax31 = Axis(f3[1,1], xlabel = "energy (eV)", ylabel = "energy bin size (eV)",
          title = "Energy variation of energy bin size")
scatterlines!(E, dE, color = :red)
display(f3)

## Plot cross-sections and back-scattering ratios.
# This will have to be done, eventually... But it is not the most important at the moment
