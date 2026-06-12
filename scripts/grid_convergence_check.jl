# Quantify the spatial discretisation error of the altitude grid for a given simulation
# configuration. Runs the same steady-state case on the default grid and on a refined grid
# (scale = 0.25, i.e. 4x more points), then reports the differences on the QO1S
# volume-excitation profile. The transport scheme converges at first order in the step
# size, so the error of any other `scale` follows by proportionality.
#
# Run from the package root. Adjust the parameters below to the case you want to check.

using AURORA
using NCDatasets

altitude_lims = [100, 400]      # (km) altitude limits of the ionosphere
θ_lims = 180:-30:0              # (°) angle-limits for the electron beams
E_max = 500                     # (eV) upper limit to the energy grid
B_angle_to_zenith = 13          # (°) angle between the B-field line and the zenith
scale_truth = 0.25              # refinement used as the reference ("truth") grid

msis_file = "test/regression/reference_results/msis_20051008-2200_70N-19E.txt"
iri_file = "test/regression/reference_results/iri_20051008-2200_70N-19E.txt"
make_flux() = InputFlux(FlatSpectrum(1e-2; E_min=E_max - 100); beams=1:2)

function run_case(scale)
    grid = AltitudeGrid(altitude_lims[1], altitude_lims[2]; scale)
    model = AuroraModel(grid, θ_lims, E_max, msis_file, iri_file, B_angle_to_zenith)
    savedir = mktempdir()
    sim = AuroraSimulation(model, make_flux(), savedir; mode=SteadyStateMode())
    run!(sim; verbose=false)
    make_volume_excitation_file(sim)
    res = load_volume_excitation(savedir)
    return vec(res.h_atm), vec(res.QO1S[:, 1])
end

# Linear interpolation of profile (x, y) at the query altitude xq
function interp1(x, y, xq)
    i = clamp(searchsortedlast(x, xq), 1, length(x) - 1)
    w = (xq - x[i]) / (x[i + 1] - x[i])
    return (1 - w) * y[i] + w * y[i + 1]
end

column(z, Q) = sum(0.5 .* (Q[1:end-1] .+ Q[2:end]) .* diff(z))

z_truth, Q_truth = run_case(scale_truth)
z_def, Q_def = run_case(1.0)

peak_truth = maximum(Q_truth)
Qi = [interp1(z_def, Q_def, zq) for zq in z_truth]

println("n_z: truth (scale=$scale_truth) = ", length(z_truth), ", default = ", length(z_def))
println("peak QO1S error:   ",
        round(100 * abs(maximum(Q_def) - peak_truth) / peak_truth, sigdigits=3), " %")
println("column QO1S error: ",
        round(100 * abs(column(z_def, Q_def) - column(z_truth, Q_truth)) / column(z_truth, Q_truth),
              sigdigits=3), " %")
# Pointwise errors blow up where the profile is steep and small (the lower flank of the
# deposition layer), so report them above increasing significance thresholds
for thr in (0.01, 0.10, 0.5)
    m = Q_truth .> thr * peak_truth
    println("max pointwise error where QO1S > ", Int(100 * thr), "% of peak: ",
            round(100 * maximum(abs.(Qi[m] .- Q_truth[m]) ./ Q_truth[m]), sigdigits=3), " %")
end
