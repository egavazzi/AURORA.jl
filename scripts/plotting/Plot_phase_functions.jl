using AURORA
using CairoMakie
CairoMakie.activate!()

altitude_lims = [100, 600];
θ_lims = 180:-10:0
E_max = 3000;

msis_file = find_msis_file(year = 2005, month = 10, day = 8, hour = 22, minute = 0,
                              lat = 70, lon = 19, height = 85:1:700);
iri_file = find_iri_file(year = 2005, month = 10, day = 8, hour = 22, minute = 0, lat = 70,
                         lon = 19, height = 85:1:700);

model = AuroraModel(altitude_lims, θ_lims, E_max, msis_file, iri_file);

## Plot phase functions
figs = plot_model(model; panels = [:phase_functions])
display(figs[:phase_functions])