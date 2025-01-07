# Function to restructure the matrix from 3D [n_mu x nz, nt, nE] to 4D [n_mu, nz, nt, nE]
function restructure_Ie_from_3D_to_4D(Ie_raw, μ_lims, h_atm, t_run, E)
    n_μ = length(μ_lims) - 1
    n_z = length(h_atm)
    n_t = length(t_run)
    n_E = length(E)
    Ie_restructured = zeros(n_μ, n_z, n_t, n_E);
    for i_E in 1:n_E
        for i_t in 1:n_t
            for i_z in 1:n_z
                for i_μ in 1:n_μ
                    Ie_restructured[i_μ, i_z, i_t, i_E] = Ie_raw[i_z + (i_μ - 1) * n_z, i_t, i_E]
                end
            end
        end
    end
    return Ie_restructured # size [n_mu, nz, nt, nE]
end


"""
    restructure_streams_of_Ie(Ie, θ_lims, new_θ_lims)

Function that merges the streams of `Ie` that are given over `θ_lims` to fit the
`new_θ_lims` of interest. It can be useful when wanting to merge some streams for plotting.

For example, if we have
```
    θ_lims = [180 160 140 120 100 90 80 60 40 20 0] # simulation
```
and we want to plot with (this should be an array of tuples, but to simplify the comparison we write it as a vector here)
```
    new_θ_lims = [180 160 120 100 90 80 60 40 20 0] # to plot
```
the function will merge the streams (160°-140°) and (140°-120°) together into a new stream with limits (160°-120°).

# Calling
`Ie_plot = restructure_streams_of_Ie(Ie, θ_lims, new_θ_lims)`

# Arguments
- `Ie`: array of electron flux with pitch-angle limits `θ_lims`. Of shape [n\\_μ, n\\_z, n\\_t, n\\_E].
- `θ_lims`: pitch-angle limits. Usually a vector or range.
- `new_θ_lims`: new pitch-angle limits. Given as an array of tuples with two rows, for example:
```
            julia> new_θ_lims = [(0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90);  # DOWN
                                 (0, 10)   (10, 30)   (30, 60)   (60, 80)   (80, 90)]  # UP
```

# Returns
- `Ie_plot`: array of electron flux with the new pitch-angle limits `new_θ_lims`. Of shape
             [n\\_μ\\_new, n\\_z, n\\_t, n\\_E], where n\\_μ\\_new is the number of streams
             in `new_θ_lims`. The first dimension of `Ie_plot` is sorted such that the
             indices go along the first row of `new_θ_lims`, and then the second row.
             In our example with `new_θ_lims` from above, that would be ``[1 2 3 4 5; 6 7 8 9 10]``.

"""
function restructure_streams_of_Ie(Ie, θ_lims, new_θ_lims)
    # Initialize the new Ie_plot that will contain the restructured streams
    n_μ_new = length(new_θ_lims)
    n_z = size(Ie, 2)
    n_t = size(Ie, 3)
    n_E = size(Ie, 4)
    Ie_plot = zeros(n_μ_new, n_z, n_t, n_E)

    # Modify the values of the down-angles so that field aligned is 180°.
    # For example new_θ_lims would go from something like
    #   2×5 Matrix{Tuple{Int64, Int64}}:
    #   (0, 10)  (10, 30)  (30, 60)  (60, 80)  (80, 90) # DOWN
    #   (0, 10)  (10, 30)  (30, 60)  (60, 80)  (80, 90) # UP
    # to
    #   2×5 Matrix{Tuple{Int64, Int64}}:
    #   (180, 170)  (170, 150)  (150, 120)  (120, 100)  (100, 90) # DOWN
    #   (0, 10)     (10, 30)    (30, 60)    (60, 80)    (80, 90)  # UP
    new_θ_lims_temp = copy(new_θ_lims)
    new_θ_lims_temp[1, :] = map(x -> 180 .- x, new_θ_lims_temp[1, :] )

    # Restructure `new_θ_lims_temp` from a 2D array to a 1D vector, by concatenating the rows.
    # Following the example from above, new_θ_lims would now be
    #   10-element Vector{Tuple{Int64, Int64}}:
    #   (180, 170)
    #   (170, 150)
    #   (150, 120)
    #   (120, 100)
    #   (100, 90)
    #   (0, 10)
    #   (10, 30)
    #   (30, 60)
    #   (60, 80)
    #   (80, 90)
    new_θ_lims_temp = vcat(eachrow(new_θ_lims_temp)...)

    # Restructure to [n_μ_new, n_z, n_t, n_E]
    # Loop over the new_θ_lims streams
    @views for i in eachindex(new_θ_lims_temp)
        # Find the indices of the streams from the simulation that should be merged in the stream new_θ_lims[i].
        idx_θ = axes(Ie, 1)[minimum(new_θ_lims_temp[i]) .<= acosd.(mu_avg(θ_lims)) .<= maximum(new_θ_lims_temp[i])]
        # Loop over these streams and add them into the right stream of Ie_plot.
        for j in idx_θ
            Ie_plot[i, :, :, :] .+= Ie[j, :, :, :]
        end
    end

    # # Extracts the values, remove the doublets and sort.
    # # Continuing the example from above, this gives us
    # #       [180 170 150 120 100 90 80 60 30 10 0]
    # new_θ_lims_temp = sort(unique(collect(Iterators.flatten(new_θ_lims_temp))); rev=true)

    return Ie_plot
end
