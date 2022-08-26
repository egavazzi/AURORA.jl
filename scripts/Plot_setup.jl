using CairoMakie
using Aurora

altitude_max = 400;
θ_lims = 180:-10:0;


h_atm, nN2, nO2, nO, ne, Te, E, dE, 
    N2_levels, O2_levels, O_levels, 
    σ_N2, σ_O2, σ_O,
    secondary_e_N2, secondary_e_O2, secondary_e_O,
    θ_lims, μ_lims, μ_center, μ_scatterings = setup(altitude_max, θ_lims);

## Plot densities
f1 = Figure()
ax11 = Axis(f1[1,1], xlabel = "(m⁻³)", ylabel = "height (km)", 
          title = "MSIS neutral density", xscale = log10, xticks = LogTicks(LinearTicks(9)))
lines!(nN2, h_atm/1e3, linewidth = 2, label = "N₂")
lines!(nO2, h_atm/1e3, linestyle = :dash, linewidth = 2, color = :red, label = "O₂")
lines!(nO, h_atm/1e3, linestyle = :dashdot, linewidth = 2, label = "O")
xlims!(1e12, 1e20)
ylims!(85, 400)
axislegend("Densities", position = :rt)
ax12 = Axis(f1[1,2], xlabel = "(m⁻³)", ylabel = "height (km)", 
        title = "electron density", xscale = log10)
lines!(ne, h_atm/1e3, linewidth = 2)
xlims!(1e8, 1e13)
ylims!(85, 400)
display(f1)

## Plot excitation thresholds
f2 = Figure()
ax21 = Axis(f2[1,1], ylabel = "Excitation energy (eV)", title = "N₂")
ax22 = Axis(f2[1,2], title = "O₂")
ax23 = Axis(f2[1,3], title = "O")
for i in eachindex(N2_levels[:,1])
lines!(ax21, [0, 1], N2_levels[i, 1].*[1, 1], linewidth = 2)
end
for i in eachindex(O2_levels[:,1])
lines!(ax22, [0, 1], O2_levels[i, 1].*[1, 1], linewidth = 2)
end
for i in eachindex(O_levels[:,1])
lines!(ax23, [0, 1], O_levels[i, 1].*[1, 1], linewidth = 2)
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
# This will have to be done, eventually... But it is not the most important atm