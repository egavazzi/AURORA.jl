using CairoMakie

phfcnE_N2, phfcnI_N2 = phase_fcn_N2(deg2rad.(0:180), E);
phfcnE_O2, phfcnI_O2 = phase_fcn_O2(deg2rad.(0:180), E);
phfcnE_O, phfcnI_O = phase_fcn_O(deg2rad.(0:180), E);

##
custom_formatter(values) = map(
		v -> "10" * Makie.UnicodeFun.to_superscript(round(Int64, v)),
		values)
##
set_theme!(
    Axis = (yticks = 0:30:180,
            xscale = log10, 
            xticks = LogTicks([1, 2, 3, 4]),
            xminorticksvisible = true,
            xminorticks = IntervalsBetween(9),
            limits = (5, 1e4, 0, 180)),
    Heatmap = (colormap = :jet, colorrange = (-6, -1))
)
##
f = Figure(resolution = (1000, 1000))
axE_N2 = Axis(f[1,1], title = "N₂ elastic collisions",
            ylabel = "scattering angle (°)")
hm = heatmap!(E, 0:180, log10.(phfcnE_N2)')
cb = Colorbar(f[1, 2], hm, tickformat = custom_formatter)

axI_N2 = Axis(f[1,3], title = "N₂ inelastic collisions")
hm = heatmap!(E, 0:180, log10.(phfcnI_N2)')
cb = Colorbar(f[1, 4], hm, tickformat = custom_formatter)

axE_O2 = Axis(f[2,1], title = "O₂ elastic collisions",
            ylabel = "scattering angle (°)")
hm = heatmap!(E, 0:180, log10.(phfcnE_O2)')
cb = Colorbar(f[2, 2], hm, tickformat = custom_formatter)

axI_O2 = Axis(f[2,3], title = "O₂ inelastic collisions")
hm = heatmap!(E, 0:180, log10.(phfcnI_O2)')
cb = Colorbar(f[2, 4], hm, tickformat = custom_formatter)

axE_O = Axis(f[3,1], title = "O elastic collisions",
            ylabel = "scattering angle (°)",
            xlabel = "Energy (eV)")
hm = heatmap!(E, 0:180, log10.(phfcnE_O)')
cb = Colorbar(f[3, 2], hm, tickformat = custom_formatter)

axI_O = Axis(f[3,3], title = "O inelastic collisions",
            xlabel = "Energy (eV)")
hm = heatmap!(E, 0:180, log10.(phfcnI_O)')
cb = Colorbar(f[3, 4], hm, tickformat = custom_formatter)
display(f)
##
set_theme!()