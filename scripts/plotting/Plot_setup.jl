using AURORA
using CairoMakie
CairoMakie.activate!(inline = true)

altitude_lims = [100, 600];     # (km) altitude limits of the ionosphere
θ_lims = 180:-10:0          # (°) angle-limits for the electron beams
E_max = 5000;               # (eV) upper limit to the energy grid

msis_file = find_msis_file(year = 2005, month = 10, day = 8, hour = 22, minute = 0,
                              lat = 70, lon = 19, height = 85:1:700);
iri_file = find_iri_file(year = 2005, month = 10, day = 8, hour = 22, minute = 0, lat = 70,
                         lon = 19, height = 85:1:700);

model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file);

## Plot all model panels
figs = plot_model(model)

## Display all figures
for (name, fig) in figs
    display(fig)
end
