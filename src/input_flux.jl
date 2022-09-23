using MAT

function Ie_top_from_file(filename, μ_center, t, E, n_loop)
    Ie_top = Array{Float64}(undef, length(μ_center), (n_loop - 1) * (length(t) - 1) + length(t), length(E))

    file = matopen(filename)
    Ie_top_raw = read(file, "Ie_total")
    close(file)

    # for constant input flux (e.g. first run), we need to resize the matrix from
    # [1, n_μ * [n_E, 1]] to [1, n_μ * [n_E, n_t]]
    if size(Ie_top[1], 2) == 1  
        for i_μ in eachindex(μ_center)
            Ie_top_raw[i_μ] = Ie_top_raw[i_μ] * ones(1, length(t) + (n_loop - 1) * (length(t) - 1)) # (e-/m²/s)
        end
    end

    # then we resize the matrix from [1, n_μ * [n_E, n_t]] to [n_μ, n_t, n_E] to be consistent with
    # the other flux matrices
    for i_μ in eachindex(μ_center)
        Ie_top[i_μ, :,  :] = Ie_top_raw[i_μ][1:length(E), :]' # (e-/m²/s)
    end

    return Ie_top
end