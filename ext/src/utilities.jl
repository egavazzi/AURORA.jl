# ======================================================================================== #
# Shared plotting helpers (used across plot_emission!, plot_ionization!, ...)
# ======================================================================================== #

"""
    _select_time_index(t, n_t, time_index) -> (i_t, label)

Return the time column index and a descriptive label for profile-mode plotting.
- If `time_index` is `nothing`, assumes steady-state and returns `(1, "Steady state")`.
- Otherwise validates bounds and returns `(time_index, "t = X.XXX s")`.
"""
function _select_time_index(t::AbstractVector, n_t::Int, time_index::Union{Int,Nothing})
    if time_index === nothing
        return 1, "Steady state"
    else
        1 <= time_index <= n_t || throw(ArgumentError(
            "time_index=$time_index is out of range [1, $n_t]"))
        return time_index, "t = $(round(t[time_index]; digits=3)) s"
    end
end

"""
    _profile_xlims_from_top(q)

Compute x-limits for log-x profile plots. The lower limit is `q[end-3] / 10`
(a few points from the top to avoid 0 values around boundary) and the upper limit is
`maximum(q) * 10`.
"""
function _profile_xlims_from_top(q::AbstractVector)
    xmin = q[end - 3] / 10
    xmax = maximum(q) * 10
    return (xmin, xmax)
end

"""
    _height_of_max(Q, h_atm)

Compute the altitude of maximum emission/ionization for each time step.
Values are set to NaN where the maximum is less than 1/10 of the global maximum.
"""
function _height_of_max(Q, h_atm)
    h_max = [h_atm[i_max[1]] for i_max in vec(findmax(Q, dims=1)[2])]
    h_max[vec(maximum(Q, dims=1)) .< maximum(Q) / 10] .= NaN
    return h_max
end
