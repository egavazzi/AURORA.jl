# Compare the QO1S volume-excitation profiles of two analysis files
# (savedir/analysis/volume_excitation.nc), e.g. before/after a numerical change or two
# runs of the same case on different altitude grids. The profiles may be on different
# altitude grids: file B is interpolated onto the altitudes of file A.
#
# Run from the package root. Set the two file paths below.

using NCDatasets

file_A = "test/regression/reference_results/SS/volume_excitation.nc"  # reference
file_B = "..."                                                        # comparison

function load_profile(file)
    NCDataset(file, "r") do ds
        vec(Array(ds["altitude"])), Array(ds["QO1S"])
    end
end

# Linear interpolation of profile (x, y) at the query altitude xq
function interp1(x, y, xq)
    i = clamp(searchsortedlast(x, xq), 1, length(x) - 1)
    w = (xq - x[i]) / (x[i + 1] - x[i])
    return (1 - w) * y[i] + w * y[i + 1]
end

column(z, Q) = sum(0.5 .* (Q[1:end-1] .+ Q[2:end]) .* diff(z))

function compare(z_A, Q_A, z_B, Q_B; label="")
    peak_A, i_peak = findmax(Q_A)
    println("-- ", label, " --")
    println("peak QO1S: A = ", peak_A, " at z = ", round(z_A[i_peak] / 1e3, digits=1),
            " km, B = ", maximum(Q_B))
    println("peak difference:   ",
            round(100 * abs(maximum(Q_B) - peak_A) / peak_A, sigdigits=3), " %")
    println("column difference: ",
            round(100 * abs(column(z_B, Q_B) - column(z_A, Q_A)) / column(z_A, Q_A),
                  sigdigits=3), " %")
    Qi = [interp1(z_B, Q_B, zq) for zq in z_A]
    # Pointwise differences blow up where the profile is steep and small (the lower flank
    # of the deposition layer), so report them above increasing significance thresholds
    for thr in (0.01, 0.10, 0.5)
        m = Q_A .> thr * peak_A
        println("max pointwise difference where QO1S > ", Int(100 * thr), "% of peak: ",
                round(100 * maximum(abs.(Qi[m] .- Q_A[m]) ./ Q_A[m]), sigdigits=3), " %")
    end
end

z_A, Q_A = load_profile(file_A)
z_B, Q_B = load_profile(file_B)
println("n_z: A = ", length(z_A), ", B = ", length(z_B),
        "  n_t: A = ", size(Q_A, 2), ", B = ", size(Q_B, 2))

if size(Q_A, 2) == 1 && size(Q_B, 2) == 1
    compare(z_A, vec(Q_A[:, 1]), z_B, vec(Q_B[:, 1]); label="steady state")
else
    # Time-dependent runs: compare the time-averaged profile and the final-time profile
    compare(z_A, vec(sum(Q_A; dims=2)) ./ size(Q_A, 2),
            z_B, vec(sum(Q_B; dims=2)) ./ size(Q_B, 2); label="time average")
    compare(z_A, vec(Q_A[:, end]), z_B, vec(Q_B[:, end]); label="final time")
end
