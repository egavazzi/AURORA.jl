# using PythonCall
import DataInterpolations

function e_O2elastic(Ep)
    s = [4.0725, 4.2883, 5.272, 5.8454, 7.1862, 8.8346, 9.7956, 10.861, 8.8346, 6.4813, 4.5155, 2.2608, 0.91126, 0.10422] .* 1e-20
    E = [0.05, 0.1, 0.6, 1, 2, 7, 10, 13.58, 30, 60, 100, 300, 1000, 1e4]

    # pyinterpolate = pyimport("scipy.interpolate")
    # cross_section = [pyinterpolate.PchipInterpolator(log.(E), log.(s))(log.(Ep[Ep .< E[end]])),
    #                  pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    # cross_section = pyconvert.(Array, cross_section)
    # cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = [
        DataInterpolations.PCHIPInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Extension)(log.(Ep[Ep .< E[end]]));
        DataInterpolations.LinearInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Linear)(log.(Ep[Ep .>= E[end]]))
        ]
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 1e-2] .= 0

    return cross_section
end

function e_O2_OO3S(Ep)
    s = [0.15583, 1.2885, 2.185, 2.9369, 3.3692, 3.3338, 2.6989, 1.7689] .* 1e-22
    E = [15.647, 18.361, 30, 50, 80, 100, 200, 400]

    # pyinterpolate = pyimport("scipy.interpolate")
    # cross_section = [pyinterpolate.PchipInterpolator(log.(E), log.(s))(log.(Ep[Ep .< E[end]])),
    #                  pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    # cross_section = pyconvert.(Array, cross_section)
    # cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = [
        DataInterpolations.PCHIPInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Extension)(log.(Ep[Ep .< E[end]]));
        DataInterpolations.LinearInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Linear)(log.(Ep[Ep .>= E[end]]))
        ]
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 15.6] .= 0

    return cross_section
end

function e_O2_9p97(Ep)
    s = [0.12616, 1.4321, 3.7052, 5.9598, 6.9829, 6.283, 3.5146, 2.185, 1.316] .* 1e-22
    E = [10.66, 20, 40, 60, 90, 100, 200, 400, 600]

    # pyinterpolate = pyimport("scipy.interpolate")
    # cross_section = [pyinterpolate.PchipInterpolator(log.(E), log.(s))(log.(Ep[Ep .< E[end]])),
    #                  pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    # cross_section = pyconvert.(Array, cross_section)
    # cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = [
        DataInterpolations.PCHIPInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Extension)(log.(Ep[Ep .< E[end]]));
        DataInterpolations.LinearInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Linear)(log.(Ep[Ep .>= E[end]]))
        ]
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 9.97] .= 0

    return cross_section
end

function e_O2_8p4(Ep)
    s = [0.68369, 9.5863, 97.91, 120.94, 120.94, 97.91, 44.339, 13.16] .* 1e-22
    E = [8.7992, 9, 10, 20, 50, 90, 120, 600]

    # pyinterpolate = pyimport("scipy.interpolate")
    # cross_section = [pyinterpolate.PchipInterpolator(log.(E), log.(s))(log.(Ep[Ep .< E[end]])),
    #                  pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    # cross_section = pyconvert.(Array, cross_section)
    # cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = [
        DataInterpolations.PCHIPInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Extension)(log.(Ep[Ep .< E[end]]));
        DataInterpolations.LinearInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Linear)(log.(Ep[Ep .>= E[end]]))
        ]
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 8.67] .= 0

    return cross_section
end

function e_O2_6(Ep)
    s = [59.598, 146.26, 235.27, 223.16, 200.79, 56.532, 45.767, 28.453, 19.66, 9.3859, 6.8369] .* 1e-23
    E = [6, 7, 8, 9, 10, 17.783, 20, 30, 40, 90, 100]

    # pyinterpolate = pyimport("scipy.interpolate")
    # cross_section = [pyinterpolate.PchipInterpolator(log.(E), log.(s))(log.(Ep[Ep .< E[end]])),
    #                  pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    # cross_section = pyconvert.(Array, cross_section)
    # cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = [
        DataInterpolations.PCHIPInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Extension)(log.(Ep[Ep .< E[end]]));
        DataInterpolations.LinearInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Linear)(log.(Ep[Ep .>= E[end]]))
        ]
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 6] .= 0

    return cross_section
end

function e_O2_4p5(Ep)
    s = [4.9802, 6.8369, 56.532, 95.863, 95.863, 69.829, 7.5986, 1.7319] .* 1e-23
    E = [4.9, 5, 6, 7, 8, 10, 20, 30]

    # pyinterpolate = pyimport("scipy.interpolate")
    # cross_section = [pyinterpolate.PchipInterpolator(log.(E), log.(s))(log.(Ep[Ep .< E[end]])),
    #                  pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    # cross_section = pyconvert.(Array, cross_section)
    # cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = [
        DataInterpolations.PCHIPInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Extension)(log.(Ep[Ep .< E[end]]));
        DataInterpolations.LinearInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Linear)(log.(Ep[Ep .>= E[end]]))
        ]
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 4.9] .= 0

    return cross_section
end

function e_O2b1Sgp(Ep)
    s = [1.9249, 2.6425, 9.3859, 15.916, 19.66, 19.66, 15.916, 9.0932, 5.535, 1.0213, 0.28755] .* 1e-23
    E = [1.9, 2, 3, 4, 5, 7, 10, 20, 30, 90, 120]

    # pyinterpolate = pyimport("scipy.interpolate")
    # cross_section = [pyinterpolate.PchipInterpolator(log.(E), log.(s))(log.(Ep[Ep .< E[end]])),
    #                  pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    # cross_section = pyconvert.(Array, cross_section)
    # cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = [
        DataInterpolations.PCHIPInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Extension)(log.(Ep[Ep .< E[end]]));
        DataInterpolations.LinearInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Linear)(log.(Ep[Ep .>= E[end]]))
        ]
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 1.9] .= 0

    return cross_section
end

function e_O2a1Dg(Ep)
    s = [6.1516, 15.916, 37.052, 56.532, 73.616, 86.254, 89.031, 75.986, 33.338, 19.66, 5.535, 2.9369] .* 1e-23
    E = [1.4678, 2, 3, 4, 5, 6, 7, 9, 20, 30, 80, 100]

    # pyinterpolate = pyimport("scipy.interpolate")
    # cross_section = [pyinterpolate.PchipInterpolator(log.(E), log.(s))(log.(Ep[Ep .< E[end]])),
    #                  pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    # cross_section = pyconvert.(Array, cross_section)
    # cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = [
        DataInterpolations.PCHIPInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Extension)(log.(Ep[Ep .< E[end]]));
        DataInterpolations.LinearInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Linear)(log.(Ep[Ep .>= E[end]]))
        ]
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 0.977] .= 0

    return cross_section
end

function e_O2vib(Ep)
    s = [1.9973, 17.464, 26.395, 23.806, 1.255, 2.5856, 11.555, 1.6247] .* 1e-22
    E = [0.24755, 0.44367, 0.6389, 0.73923, 1.4251, 6.1282, 9.086, 14.7]

    # pyinterpolate = pyimport("scipy.interpolate")
    # cross_section = [pyinterpolate.PchipInterpolator(log.(E), log.(s))(log.(Ep[Ep .< E[end]])),
    #                  pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    # cross_section = pyconvert.(Array, cross_section)
    # cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = [
        DataInterpolations.PCHIPInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Extension)(log.(Ep[Ep .< E[end]]));
        DataInterpolations.LinearInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Linear)(log.(Ep[Ep .>= E[end]]))
        ]
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 0.24] .= 0

    return cross_section
end

function e_O2ionx2pg(Ep)
    # e_O2ionX2Pg - electron ionisation cross section (m^2) to the
    # ground-state of O2+
    cross_section = similar(Ep)
    for iE in eachindex(Ep)
        if Ep[iE] > 12.1 && Ep[iE] < 300
            cross_section[iE] = (1-12.1/Ep[iE])*exp(-92.26781+44.59405*log(Ep[iE])-13.45450*log(Ep[iE])^2+1.809351*log(Ep[iE])^3-0.09207181*log(Ep[iE])^4)
        elseif Ep[iE] >= 300
            cross_section[iE] = 8.771133e-15*log(0.024825*Ep[iE])/Ep[iE]
        else
            cross_section[iE] = 0
        end
    end
    return cross_section ./ 1e4
end

function e_O2iona4pu(Ep)
    # e_O2iona4Pu - electron ionisation cross section (m^2) to the
    # first excited state of O2+
    # Ep electron energy (eV)
    cross_section = similar(Ep)
    for iE in eachindex(Ep)
        if Ep[iE] > 16.1 && Ep[iE] < 300
            cross_section[iE] = (1-16.1/Ep[iE])*exp(-49.02795+6.907601*log(Ep[iE])-1.325630*log(Ep[iE])^2+0.08157475*log(Ep[iE])^3-0.0002650063*log(Ep[iE])^4)
        elseif Ep[iE] >= 300
            cross_section[iE] = 4.911588e-15*log(0.024825*Ep[iE])/Ep[iE]
        else
            cross_section[iE] = 0
        end
    end
    return cross_section ./ 1e4
end

function e_O2ion16p9(Ep)
    # e_O2ion16p9 - electron ionisation cross section (m^2) to the
    # excited states A2Pu, B2Sg-, 2Pu, c4Sg- of O2+
    # Ep electron energy (eV)
    cross_section = similar(Ep)
    for iE in eachindex(Ep)
        if Ep[iE] > 16.1 && Ep[iE] < 300
            cross_section[iE] = (1-16.1/Ep[iE])*exp(-49.02795+6.907601*log(Ep[iE])-1.325630*log(Ep[iE])^2+0.08157475*log(Ep[iE])^3-0.0002650063*log(Ep[iE])^4)
        elseif Ep[iE] >= 300
            cross_section[iE] = 4.911588e-15*log(0.024825*Ep[iE])/Ep[iE]
        else
            cross_section[iE] = 0
        end
    end
    return cross_section ./ 1e4
end

function e_O2ionb4sgm(Ep)
    # e_O2ionb4sgm - electron ionisation cross section (m^2) to the
    # excited states O_2^+(b^4\Sigma_g^-)
    # Ep electron energy (eV)
    cross_section = similar(Ep)
    for iE in eachindex(Ep)
        if Ep[iE] > 18.2 && Ep[iE] < 250
            cross_section[iE] = (1-18.2/Ep[iE])*exp(-26.45362-12.67161*log(Ep[iE])+4.799453*log(Ep[iE])^2-0.7680145*log(Ep[iE])^3+0.04377051*log(Ep[iE])^4)
        elseif Ep[iE] >= 250
            cross_section[iE] = 1.980794e-15*log(0.024825*Ep[iE])/Ep[iE]
        else
            cross_section[iE] = 0
        end
    end
    return cross_section ./ 1e4
end

function e_O2dion(Ep)
    # e_O2dion - dissociative ionization cross section (m^2)
    # Ep electron energy (eV)
    cross_section = similar(Ep)
    for iE in eachindex(Ep)
        if Ep[iE] > 18.9 && Ep[iE] < 200
            cross_section[iE] = (1-18.9/Ep[iE])*exp(-112.0394+59.82636*log(Ep[iE])-17.95648*log(Ep[iE])^2+2.411165*log(Ep[iE])^3-0.1228601*log(Ep[iE])^4)
        elseif Ep[iE] >= 200
            cross_section[iE] = 6.078054e-15*log(0.030992*Ep[iE])/Ep[iE]
        else
            cross_section[iE] = 0
        end
    end
    return cross_section ./ 1e4
end

function e_O2ddion(Ep)
    # e_O2ddion - dissociative-double ionization cross section (m^2)
    # Ep electron energy (eV)
    cross_section = similar(Ep)
    for iE in eachindex(Ep)
        if Ep[iE] > 32.51 && Ep[iE] < 300
            cross_section[iE] = (1-32.51/Ep[iE])*exp(-228.0913+141.5861*log(Ep[iE])-40.00023*log(Ep[iE])^2+5.079899*log(Ep[iE])^3-0.2454952*log(Ep[iE])^4)
        elseif Ep[iE] >= 200
            cross_section[iE] = 8.141987e-16*log(0.12511*Ep[iE])/Ep[iE]
        else
            cross_section[iE] = 0
        end
    end
    return cross_section ./ 1e4
end
