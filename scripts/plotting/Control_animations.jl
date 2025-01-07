#=
Functions to animate Ie(z, t, E).

So far we have
- animate_IeztE_3Dzoft: Ie as a heatmap over height and energy, animation in time.


=#
using AURORA
using Printf
using GLMakie
GLMakie.activate!()


## Calling the animate function

directory_to_process = "Visions2/Alfven_475s"
angles_to_plot = [(0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90);  # DOWN
                  (0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90)]  # UP
color_limits = (1e5, 1e9)
animate_IeztE_3Dzoft(directory_to_process, angles_to_plot, color_limits; plot_Ietop = true)






## Testing

# First load Ie
using MAT
filename = "data/Visions2/Alfven_475s/IeFlickering-01.mat"

data = matread(filename);
Ie_raw = data["Ie_ztE"]; # size of [n_mu x nz, nt, nE]
μ_scatterings = data["mu_scatterings"]
μ_lims = data["mu_lims"]
t_run = data["t_run"]
E = data["E"]
h_atm = data["h_atm"]
dE = diff(E); dE = [dE; dE[end]];
θ_lims = acosd.(μ_lims)

# Restructure and plot
angles_to_plot = [(0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90);  # DOWN
                 (0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90)]  # UP
Ie = AURORA.restructure_Ie_from_3D_to_4D(Ie_raw, μ_lims, h_atm, t_run, E); # size [n_mu, nz, nt, nE]
Ie_plot = AURORA.restructure_streams_of_Ie(Ie, θ_lims, angles_to_plot);

Ie_timeslice = Observable(Ie_plot[:, :, 1, :])
time = Observable("$(t_run[1]) s")
f = AURORA.make_IeztE_3Dzoft_plot(Ie_timeslice, time, h_atm, E, angles_to_plot, (1e5, 1e9), t_top, data_Ietop)
display(f)


## Test the animation
@time for i_t in 1:length(t_run)
    Ie_timeslice[] .= Ie_plot[:, :, i_t, :]
    time[] = "$(t_run[i_t]) s"
    notify(Ie_timeslice)
    sleep(0.01)
end

## Testing recording the animation to a mp4
# with glmakie(~12s for n_t=51)
@time record(f, "time_animation1.mp4", 1:10; framerate=10, backend=GLMakie) do i_t
    Ie_timeslice[] .= Ie_plot[:, :, i_t, :]
    time[] = "$(t_run[i_t]) s"
    notify(Ie_timeslice)
end
# with cairomakie (~110s for n_t=51)
@time record(f, "time_animation2.mp4", 1:10; framerate=10, backend=CairoMakie) do i_t
    Ie_timeslice[] .= Ie_plot[:, :, i_t, :]
    time[] = "$(t_run[i_t]) s"
    notify(Ie_timeslice)
end
