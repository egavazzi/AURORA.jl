using Makie

#=
This is the function that creates the plot from the data observable. Then when we update the
data, the plot is automatically updated.

It should take Ie as a function of pitch-angle, height and energy as input. But NOT as a
function of time. This because we want to update the time OUTSIDE of the function.

It is also possible to plot the precipitating Ie at the top of the ionosphere by giving an
optional NamedTuple `Ietop_struct` as input. By default it is set to
```
Ietop_struct = (bool = false, t_top = nothing, data_Ietop = nothing)
```
which won't plot anything.
To plot something, you need to give a `Ietop_struct` with
`Ie_top_struct.bool = true` as well as some values for `Ie_top_struct.t_top` and
`Ie_top_struct.data_Ietop`.
=#
function make_Ie_in_time_plot(Ie_timeslice::Observable{Array{Float64, 3}},
                                time::Observable{String}, h_atm, E, angles_to_plot, colorrange,
                                Ietop_struct = (bool = false, t_top = nothing, data_Ietop = nothing,
                                Ietop_angle_cone = nothing))

    # Slice the input Ie into its different pitch-angle components
    Ie_streams = Array{Observable}(nothing, length(angles_to_plot))
    for i in eachindex(angles_to_plot)
        Ie_streams[i] = @lift($Ie_timeslice[i, :, :]')
    end

    fig = Figure(size = (1500, 1000), fontsize=20)
    n_row = size(angles_to_plot, 1)
    n_col = size(angles_to_plot, 2)
    ga = fig[1:4, 1:n_col] = GridLayout()

    # Helper: check if all panels below row i in column j are nothing
    all_below_empty(i, j) = all(isnothing(angles_to_plot[k, j]) for k in (i+1):n_row)
    # Helper: check if all panels to the left of column j in row i are nothing
    all_left_empty(i, j) = all(isnothing(angles_to_plot[i, k]) for k in 1:(j-1))

    for i in axes(angles_to_plot, 1)
        for j in axes(angles_to_plot, 2)
            idx = (i - 1) * n_col  + j # goes along first row, then second row

            # Skip empty panels (nothing entries)
            isnothing(angles_to_plot[i, j]) && continue

            ax = Axis(ga[i, j], xscale = log10, xminorticks = IntervalsBetween(9),
                      xminorticksvisible = true, yticks = 100:100:600)
            heatmap!(E, h_atm / 1e3, Ie_streams[idx]; colorrange, colorscale = log10, colormap = :inferno)

            # Generate title based on angle values (>90° = DOWN, <90° = UP)
            θ1, θ2 = angles_to_plot[i, j]
            avg_angle = (θ1 + θ2) / 2
            direction = if θ1 < 90 < θ2 || θ2 < 90 < θ1
                ""
            elseif avg_angle >= 90
                "DOWN"
            else
                "UP"
            end
            ax.title = "$(Int(θ1))° - $(Int(θ2))° $direction"

            # Show x-axis labels if on bottom row OR all panels below are empty
            if i == n_row || all_below_empty(i, j)
                ax.xlabel = "Energy (eV)"
            else
                ax.xticklabelsvisible = false
            end
            # Show y-axis labels if on first column OR all panels to the left are empty
            if j == 1 || all_left_empty(i, j)
                ax.ylabel = "Altitude (km)"
            else
                ax.yticklabelsvisible = false
            end
        end
    end
    Colorbar(fig[:, end + 1]; limits = colorrange, scale = log10, label = "Ie (#e⁻/m²/s/eV/ster)", colormap = :inferno)

    # Plot Ie precipitating at the top (#TODO: this can probably be done in a better way than with a struct, or?)
    if Ietop_struct.bool
        plot_hposition = 1:floor(Int, n_col / 2)
        gb = fig[0, plot_hposition] = GridLayout()
        ax_Ietop = Axis(gb[1, 1], yscale = log10, ylabel = " Energy (eV)", xlabel = "t (s)",
                        title = "Incoming energy flux",
                        yminorticksvisible = true, yminorticks = IntervalsBetween(9),
                        xticklabelsvisible = true, xminorticksvisible = true,
                        xticksmirrored = true, yticksmirrored = true,
                        limits = ((0, 1), nothing))
        hm = heatmap!(Ietop_struct.t_top, E, Ietop_struct.data_Ietop; colormap = :inferno,
                      colorscale = log10,
                      colorrange = (1e6, maximum(Ietop_struct.data_Ietop)))
        time_float64 = @lift(parse(Float64, $time[1:end-1]))
        vlines!(time_float64, linewidth = 3)
        angles = Ietop_struct.Ietop_angle_cone
        Colorbar(gb[1, 2], hm; label = "IeE ($(180 - angles[2])°-$(180 - angles[1])°)\n(eV/m²/s/eV/ster)")

        time_hposition = ceil(Int, n_col / 2):n_col
        Label(fig[0, time_hposition], time; tellwidth = false, tellheight = false, fontsize=20)
    else
        Label(fig[0, :], time; tellwidth = false, fontsize=20)
    end

    return fig
end
