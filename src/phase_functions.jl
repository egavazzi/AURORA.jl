
function build_DCSN2()
    # Pre-allocate coefficients
    CoeffA1 = [-7.55013, 6.4651, -4.61116, 1.51264, -0.24804, 0.020829, -0.7241e-3]
    CoeffA2 = [-121.79, 161.31, -88.9543, 25.0846, -3.79789, 0.294225, -0.91716e-2]
    CoeffA3 = [-24.4482, 11.6556, -0.109564, -1.50976, 0.479943, -0.0597, 0.263263e-2]
    CoeffEt1 = [-9.06218, 11.1088, -6.1393, 1.43937, -0.149179, 0.570119e-2, 0.0]
    CoeffB1 = [3.481694, -0.47699, 0.0]
    CoeffB2 = [-42.1667, 13.79589, -1.123838]
    CoeffB3 = [-35.49183, 12.7418, -1.247576]
    CoeffEt2 = [1.590443, -0.207755, 0.0]

    # Pre-allocate working arrays
    Elog = ones(7)
    B = zeros(3)

    function DCSN2(θ, Energy)
        Elog[2] = log(Energy)
        for i in 3:7
            Elog[i] = Elog[i-1] * Elog[2]
        end

        if Energy <= 500.0
            B[1] = exp(Elog[1] * CoeffA1[1] + Elog[2] * CoeffA1[2] + Elog[3] * CoeffA1[3] +
                       Elog[4] * CoeffA1[4] + Elog[5] * CoeffA1[5] + Elog[6] * CoeffA1[6] +
                       Elog[7] * CoeffA1[7])
            Etta = exp(Elog[1] * CoeffEt1[1] + Elog[2] * CoeffEt1[2] +
                       Elog[3] * CoeffEt1[3] + Elog[4] * CoeffEt1[4] +
                       Elog[5] * CoeffEt1[5] + Elog[6] * CoeffEt1[6] +
                       Elog[7] * CoeffEt1[7])
        else
            B[1] = 2.0 - exp(Elog[1] * CoeffB1[1] + Elog[2] * CoeffB1[2] + Elog[3] * CoeffB1[3])
            Etta = 1.6038 - exp(Elog[1] * CoeffEt2[1] + Elog[2] * CoeffEt2[2])
        end

        if Energy <= 250.0
            B[2] = exp(Elog[1] * CoeffA2[1] + Elog[2] * CoeffA2[2] + Elog[3] * CoeffA2[3] +
                       Elog[4] * CoeffA2[4] + Elog[5] * CoeffA2[5] + Elog[6] * CoeffA2[6] +
                       Elog[7] * CoeffA2[7])
        else
            B[2] = exp(Elog[1] * CoeffB2[1] + Elog[2] * CoeffB2[2] + Elog[3] * CoeffB2[3])
        end

        if Energy <= 100.0
            B[3] = exp(Elog[1] * CoeffA3[1] + Elog[2] * CoeffA3[2] + Elog[3] * CoeffA3[3] +
                       Elog[4] * CoeffA3[4] + Elog[5] * CoeffA3[5] + Elog[6] * CoeffA3[6] +
                       Elog[7] * CoeffA3[7])
        else
            B[3] = exp(Elog[1] * CoeffB3[1] + Elog[2] * CoeffB3[2] + Elog[3] * CoeffB3[3])
        end

        T = Energy / 510879.0;
        Etta = Etta * 6.221e-5 / T / (T + 2);
        x = (1 .- cos.(θ) .+ 2 * Etta);

        DCS = B[1] ./ x.^2 + B[2] ./ x + B[3] .* x.^4

        return DCS
    end

    return DCSN2
end

function build_DCSO2()
    # Pre-allocate coefficients
    CoeffA1 = [-59.543, 79.8, -49.7886, 16.3775, -2.90625, 0.26045, -0.9109e-2]
    CoeffA2 = [-20.3557, 33.1387, -30.4416, 14.6647, -3.7744, 0.489333, -0.024934]
    CoeffA3 = [-15.861, 4.7481, -0.432, 0.0, 0.0, 0.0, 0.0]
    CoeffEt1 = [0.02369, 0.5232, -0.585, 0.240122, -0.041798, 0.2622e-2, 0.0]
    CoeffB1 = [3.381, -0.4622, 0.0]
    CoeffB2 = [-121.6, 41.0, -3.447]
    CoeffB3 = [-86.19, 33.04, -3.259]
    CoeffEt2 = [2.308, -0.3645, 0.0]

    # Pre-allocate working arrays
    Elog = ones(7)
    B = zeros(3)

    function DCSO2(θ, Energy)
        Elog[2] = log(Energy)
        for i in 3:7
            Elog[i] = Elog[i-1] * Elog[2]
        end

        if Energy <= 500.0
            B[1] = exp(Elog[1] * CoeffA1[1] + Elog[2] * CoeffA1[2] + Elog[3] * CoeffA1[3] +
                       Elog[4] * CoeffA1[4] + Elog[5] * CoeffA1[5] + Elog[6] * CoeffA1[6] +
                       Elog[7] * CoeffA1[7])
            Etta = (Elog[1] * CoeffEt1[1] + Elog[2] * CoeffEt1[2] + Elog[3] * CoeffEt1[3] +
                    Elog[4] * CoeffEt1[4] + Elog[5] * CoeffEt1[5] + Elog[6] * CoeffEt1[6])
        else
            B[1] = 2.0 - exp(Elog[1] * CoeffB1[1] + Elog[2] * CoeffB1[2])
            Etta = 1.3141 - exp(Elog[1] * CoeffEt2[1] + Elog[2] * CoeffEt2[2])
        end

        if Energy <= 250.0
            B[2] = exp(Elog[1] * CoeffA2[1] + Elog[2] * CoeffA2[2] + Elog[3] * CoeffA2[3] +
                       Elog[4] * CoeffA2[4] + Elog[5] * CoeffA2[5] + Elog[6] * CoeffA2[6] +
                       Elog[7] * CoeffA2[7])
        else
            B[2] = exp(Elog[1] * CoeffB2[1] + Elog[2] * CoeffB2[2] + Elog[3] * CoeffB2[3])
        end

        if Energy <= 100.0
            B[3] = exp(Elog[1] * CoeffA3[1] + Elog[2] * CoeffA3[2] + Elog[3] * CoeffA3[3])
        else
            B[3] = exp(Elog[1] * CoeffB3[1] + Elog[2] * CoeffB3[2] + Elog[3] * CoeffB3[3])
        end

        T = Energy / 510879.0
        Etta = Etta * 6.8e-5 / T / (T + 2)
        x = (1 .- cos.(θ) .+ 2 * Etta)

        DCS = B[1] ./ x.^2 + B[2] ./ x + B[3] .* x.^4

        return DCS
    end

    return DCSO2
end

function build_DCSO()
    # Pre-allocate coefficients
    CoeffA1 = [-12.61, 7.773, -3.922, 0.825, -0.05546, 0.0, 0.0]
    CoeffA2 = [-20.30, 16.37, -6.604, 1.218, -0.08263, 0.0, 0.0]
    CoeffA3 = [-16.22, 4.452, -0.3577, 0.0, 0.0, 0.0, 0.0]
    CoeffEt1 = [-4.008, 1.017, -0.06066, 0.0, 0.0, 0.0, 0.0]
    CoeffB1 = [11.44, -1.986, 0.0]
    CoeffB3 = [-19.04, 5.805, -0.5186]
    CoeffEt2 = [3.291, -0.7007, 0.0]

    # Pre-allocate working arrays
    Elog = ones(7)
    B = zeros(3)

    function DCSO(θ, Energy)
        Elog[2] = log(Energy)
        for i in 3:7
            Elog[i] = Elog[i-1] * Elog[2]
        end

        if Energy ≤ 500.0
            B[1] = exp(Elog[1] * CoeffA1[1] + Elog[2] * CoeffA1[2] + Elog[3] * CoeffA1[3] +
                       Elog[4] * CoeffA1[4] + Elog[5] * CoeffA1[5])
            Etta = exp(Elog[1] * CoeffEt1[1] + Elog[2] * CoeffEt1[2] +
                       Elog[3] * CoeffEt1[3])
        else
            B[1] = 1.0 - exp(Elog[1] * CoeffB1[1] + Elog[2] * CoeffB1[2])
            Etta = 1.3198 - exp(Elog[1] * CoeffEt2[1] + Elog[2] * CoeffEt2[2])
        end

        B[2] = exp(Elog[1] * CoeffA2[1] + Elog[2] * CoeffA2[2] + Elog[3] * CoeffA2[3] +
                   Elog[4] * CoeffA2[4] + Elog[5] * CoeffA2[5])

        if Energy ≤ 100.0
            B[3] = exp(Elog[1] * CoeffA3[1] + Elog[2] * CoeffA3[2] + Elog[3] * CoeffA3[3] +
                       Elog[4] * CoeffA3[4] + Elog[5] * CoeffA3[5])
        else
            B[3] = exp(Elog[1] * CoeffB3[1] + Elog[2] * CoeffB3[2] + Elog[3] * CoeffB3[3])
        end

        T = Energy / 510879.0
        Etta = Etta * 6.8e-5 / T / (T + 2)

        if Energy ≤ 500
            x = 4.4549 - 0.003114 * Energy - 0.4663 * Elog[2]
            x = exp(x)
            x2 = x * x
            Funθ = θ ./ (x2 - 8100) .* (x2 .- 8100 .+ x*180 .* (1 .- θ ./ 3.14159265))
        else
            Funθ = θ
        end
        x = (1 .- cos.(Funθ) .+ 2 * Etta)

        DCS = B[1] ./ x.^2 + B[2] ./ x.^3 + B[3] .* x.^4

        return DCS
    end

    return DCSO
end

const DCSN2 = build_DCSN2()
const DCSO2 = build_DCSO2()
const DCSO = build_DCSO()


## ----------------------------------------------------- ##

function phase_fcn_N2(θ, Energy)
    DCS = zeros(length(θ), length(Energy))
    for iE = length(Energy):-1:1
        for i_θ = length(θ):-1:1
            DCS[i_θ, iE] = DCSN2(θ[i_θ], Energy[iE]);
        end
    end
    i = findfirst(≥(4.796), Energy);
    if i !== nothing
        for iE in 1:(i - 1)
            DCS[:, iE] = DCS[:, i]
        end
    end
    phfcnE = DCS ./ repeat(sum(DCS, dims=1), length(θ), 1);

    for iE = length(Energy):-1:1
        for i_θ = length(θ):-1:1
            DCS[i_θ, iE] = DCSN2(θ[i_θ], Energy[iE]);
        end
    end
    phfcnI = DCS ./ repeat(sum(DCS, dims=1), length(θ), 1);

    return phfcnE, phfcnI
end

function phase_fcn_O2(θ, Energy)
    DCS = zeros(length(θ), length(Energy))
    for iE = length(Energy):-1:1
        for i_θ = length(θ):-1:1
            DCS[i_θ, iE] = DCSO2(θ[i_θ], Energy[iE]);
        end
    end
    i = findfirst(≥(4.15), Energy);
    if i !== nothing
        for iE in 1:(i-1)
            DCS[:, iE] = DCS[:, i];
        end
    end
    phfcnE = DCS ./ repeat(sum(DCS, dims=1), length(θ), 1);

    for iE = length(Energy):-1:1
        for i_θ = length(θ):-1:1
            DCS[i_θ, iE] = DCSO2(θ[i_θ], Energy[iE]);
        end
    end
    phfcnI = DCS ./ repeat(sum(DCS, dims=1), length(θ), 1);

    return phfcnE, phfcnI
end

function phase_fcn_O(θ, Energy)
    DCS = zeros(length(θ), length(Energy))
    for iE = length(Energy):-1:1
        for i_θ = length(θ):-1:1
            DCS[i_θ, iE] = DCSO(θ[i_θ], Energy[iE]);
        end
    end
    i = findfirst(≥(3.0), Energy);
    if i !== nothing
        for iE in 1:(i-1)
            DCS[:, iE] = DCS[:, i];
        end
    end
    phfcnE = DCS ./ repeat(sum(DCS, dims=1), length(θ), 1);

    for iE = length(Energy):-1:1
        if Energy[iE] <= 150
            Etmp = Energy[iE] * exp(0.025 * (Energy[iE] - 5.0));
        else
            Etmp = 50 * Energy[iE]
        end

        for i_θ = length(θ):-1:1
            DCS[i_θ, iE] = DCSO(θ[i_θ], Etmp);
        end
    end
    phfcnI = DCS ./ repeat(sum(DCS, dims=1), length(θ), 1);

    return phfcnE, phfcnI
end

function convert_phase_fcn_to_3D!(phase_fcn, θ)
    # The measurements of scattering probabilities that make up the phase function matrices were done
    # in a plane (2D). Problem with that is that the electrons scatter in the 3 dimensions, so only
    # a fraction of them are measured. This fraction depends on the scattering angle, as the e- will
    # scatter on a "ring" (think about a slice of a sphere) of area 2π*sin(θ)*dθ. This means that the
    # probability of scattering will be underestimated the more we approach angles around 90°. This
    # function is here to correct that.
    phase_fcn = phase_fcn .* abs.(sin.(θ));       # we don't need the factors 2π and dθ as they are the same for all theta bins.
    phase_fcn = phase_fcn ./ sum(phase_fcn);    # so that sum of probabilities = 1
end

function convert_phase_fcn_to_3D(phase_fcn, θ)
    # The measurements of scattering probabilities that make up the phase function matrices were done
    # in a plane (2D). Problem with that is that the electrons scatter in the 3 dimensions, so only
    # a fraction of them are measured. This fraction depends on the scattering angle, as the e- will
    # scatter on a "ring" (think about a slice of a sphere) of area 2π*sin(θ)*dθ. This means that the
    # probability of scattering will be underestimated the more we approach angles around 90°. This
    # function is here to correct that.
    phase_fcn = phase_fcn .* abs.(sin.(θ));       # we don't need the factors 2π and dθ as they are the same for all theta bins.
    phase_fcn = phase_fcn ./ sum(phase_fcn);    # so that sum of probabilities = 1
    return phase_fcn
end
