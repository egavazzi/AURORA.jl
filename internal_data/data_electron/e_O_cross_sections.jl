using DataInterpolations: LinearInterpolation, PCHIPInterpolation, ExtrapolationType

function e_Oelastic(Ep)
    # e_Oelastic - elastic electron collision cross section (m^2)
    # Ep electron energy (eV)
    cross_section = similar(Ep)
    for iE in eachindex(Ep)
        if Ep[iE] > 0.2 && Ep[iE] < 10000
            cross_section[iE] = 2.8e-17*exp(2.8+0.39*log(Ep[iE])-0.078*log(Ep[iE])^2-0.0102*log(Ep[iE])^3+0.00095*log(Ep[iE])^4)
        elseif Ep[iE] >= 10000
            cross_section[iE] = 2.8e-17*2572.3/Ep[iE]
        else
            cross_section[iE] = 0
        end
    end
    return cross_section ./ 1e4
end

function e_O1D(Ep)
    # e_O1D - electron excitation cross section (m^2)
    # Ep electron energy (eV)
    cross_section = similar(Ep)
    for iE in eachindex(Ep)
        if Ep[iE] > 1.967 && Ep[iE] < 6.867
            cross_section[iE] = (1-1.9/Ep[iE])*exp(-38.0685-0.2992*Ep[iE]+0.20375*Ep[iE]^2-0.0211739*Ep[iE]^3)
        elseif Ep[iE] > 6.867 && Ep[iE] < 30
            cross_section[iE] = exp(-34.081-0.912397*Ep[iE]+7.185417e-2*Ep[iE]^2-2.48398e-3*Ep[iE]^3+3.00574e-5*Ep[iE]^4)
        elseif Ep[iE] >= 30
            cross_section[iE] = 1.881e-13/Ep[iE]^3
        else
            cross_section[iE] = 0
        end
    end
    return cross_section ./ 1e4
end

function e_O1S(Ep)
    # e_O1S - electron excitation cross section (m^2)
    # Ep electron energy (eV)
    cross_section = similar(Ep)
    for iE in eachindex(Ep)
        if Ep[iE] > 4.17 && Ep[iE] < 30
            cross_section[iE] = (1-4.17/Ep[iE])*exp(-77.00326+61.07961*log(Ep[iE])-37.16693*log(Ep[iE])^2+10.03347*log(Ep[iE])^3-1.021318*log(Ep[iE])^4)
        elseif Ep[iE] >= 30
            cross_section[iE] = 3.24e-14/Ep[iE]^3
        else
            cross_section[iE] = 0
        end
    end
    return cross_section ./ 1e4
end

function e_O3s5S0(Ep)
    # e_O3s5S0 - electron excitation cross section (m^2)
    # Ep electron energy (eV)
    cross_section = similar(Ep)
    for iE in eachindex(Ep)
        if Ep[iE] > 9.14 && Ep[iE] < 23
            cross_section[iE] = (1-9.14/Ep[iE])*exp(-280.6036+256.7227*log(Ep[iE])-89.59541*log(Ep[iE])^2+10.22137*log(Ep[iE])^3)
        elseif Ep[iE] >= 23
            cross_section[iE] = 7.6e-15/Ep[iE]^3
        else
            cross_section[iE] = 0
        end
    end
    return cross_section ./ 1e4
end

function e_O3s3S0(Ep)
    # e_O3s3S0 - electron excitation cross section (m^2)
    # Ep electron energy (eV)
    s = [4, 9.83, 10.67, 11.17, 8.67, 7.9, 6.33, 5.67, 4.27] .* 1e-22
    E = [9.6, 13.33, 16.667, 20, 30, 50, 100, 150, 200]

    # pyinterpolate = pyimport("scipy.interpolate")
    # cross_section = [pyinterpolate.PchipInterpolator(log.(E), log.(s))(log.(Ep[Ep .< E[end]])),
    #                  pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    # cross_section = pyconvert.(Array, cross_section)
    # cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = [
        PCHIPInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Extension)(log.(Ep[Ep .< E[end]]));
        LinearInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Linear)(log.(Ep[Ep .>= E[end]]))
        ]
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 9.521] .= 0

    return cross_section
end

function e_O3p5P(Ep)
    # e_O3p5P - electron excitation cross section (m^2)
    # Ep electron energy (eV)
    cross_section = similar(Ep)
    for iE in eachindex(Ep)
        if Ep[iE] > 10.73 && Ep[iE] < 31.614
            cross_section[iE] = (1 - 10.73 / Ep[iE]) * exp(678.5487 - 931.9045 * log(Ep[iE]) + 452.4033 * log(Ep[iE])^2 - 97.15600 * log(Ep[iE])^3 + 7.763117 * log(Ep[iE])^4)
        elseif Ep[iE] > 31.614 && Ep[iE] < 200
            cross_section[iE] = exp(-54.32016 + 19.35011 * log(Ep[iE]) - 8.191310 * log(Ep[iE])^2 + 1.323729 * log(Ep[iE])^3 - 0.07969778 * log(Ep[iE])^4)
        elseif Ep[iE] >= 200
            cross_section[iE] = 1.593536e-14 / Ep[iE]^3
        else
            cross_section[iE] = 0
        end
    end
    return cross_section ./ 1e4
end

function e_O3sp3D0(Ep)
    # e_O3sp3D0 - electron excitation cross section (m^2)
    # Ep electron energy (eV)
    s = [3, 5.5, 5, 5.8, 4.5, 4, 2.5] .* 1e-22
    E = [12.6, 20, 30, 50, 100, 150, 200]

    # pyinterpolate = pyimport("scipy.interpolate")
    # cross_section = [pyinterpolate.PchipInterpolator(log.(E), log.(s))(log.(Ep[Ep .< E[end]])),
    #                  pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    # cross_section = pyconvert.(Array, cross_section)
    # cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = [
        PCHIPInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Extension)(log.(Ep[Ep .< E[end]]));
        LinearInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Linear)(log.(Ep[Ep .>= E[end]]))
        ]
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 12.54] .= 0

    return cross_section
end

function e_O3p3P(Ep)
    # e_O3p3P - electron excitation cross section (m^2)
    # Ep electron energy (eV)
    s = [5.1, 7.8, 4, 2.9, 1.1] .* 1e-22
    E = [13.5, 20, 30, 50, 100]

    # pyinterpolate = pyimport("scipy.interpolate")
    # cross_section = pyinterpolate.PchipInterpolator(log.(E), log.(s))(log.(Ep[Ep .<= 60]))
    # cross_section = pyconvert(Array, cross_section)
    cross_section = PCHIPInterpolation(log.(s), log.(E); extrapolation = ExtrapolationType.Extension)(log.(Ep[Ep .<= 60]))
    cross_section = exp.(cross_section)

    if any(Ep .> 60)
        cross_section = [cross_section; e_O3p3P([59.9]) .* 60 ./ log(60) .* log.(Ep[Ep .> 60]) ./ Ep[Ep .> 60]]
    end

    cross_section[.!isfinite.(cross_section)] .= 0
    cross_section[Ep .< 10.99] .= 0

    return cross_section
end

function e_Oion4S0(Ep)
    # e_Oion4S0 - O electron ionization cross section to O^+(4S_0) (m^2)
    # Ep - electron energy (eV)
    # Data source: Tima Sergienko, private communication
    cross_section = similar(Ep)

    for iE in eachindex(Ep)
        if Ep[iE] > 13.6 && Ep[iE] <= 250
            cross_section[iE] = 0.35*(1-13.6/Ep[iE])*exp(-38.13225-1.957729*log(Ep[iE])+1.526543*(log(Ep[iE]))^2-0.3056663*(log(Ep[iE]))^3+0.01849928*(log(Ep[iE]))^4)
        elseif Ep[iE] > 250
            cross_section[iE] = 4.760656e-15*log(0.032*Ep[iE])/Ep[iE]
        else
            cross_section[iE] = 0
        end
    end

    return cross_section ./ 1e4
end

function e_Oion2D0(Ep)
    # e_Oion2D0 - O electron ionization cross section to O^+(2P_0) (m^2)
    # Ep - electron energy (eV)
    # Data source: Tima Sergienko, private communication
    cross_section = similar(Ep)

    for iE in eachindex(Ep)
        if Ep[iE] > 16.9 && Ep[iE] <= 250
            cross_section[iE] = (1-16.9/Ep[iE])*exp(-38.80617-1.369768*log(Ep[iE])+1.103263*(log(Ep[iE]))^2-0.2226944*(log(Ep[iE]))^3+0.01331574*(log(Ep[iE]))^4)
        elseif Ep[iE] > 250
            cross_section[iE] = 4.109687e-15*log(0.032*Ep[iE])/Ep[iE]
        else
            cross_section[iE] = 0
        end
    end

    return cross_section ./ 1e4
end

function e_Oion2P0(Ep)
    # e_Oion2P0 - O electron ionization cross section to O^+(2P_0) (m^2)
    # Ep - electron energy (eV)
    # Data source: Tima Sergienko, private communication
    cross_section = similar(Ep)

    for iE in eachindex(Ep)
        if Ep[iE] > 18.6 && Ep[iE] <= 250
            cross_section[iE] = (1-18.6/Ep[iE])*exp(-34.62552-4.788984*log(Ep[iE])+2.031120*(log(Ep[iE]))^2-0.3347122*(log(Ep[iE]))^3+0.01838017*(log(Ep[iE]))^4)
        elseif Ep[iE] > 250
            cross_section[iE] = 2.330830e-15*log(0.032*Ep[iE])/Ep[iE]
        else
            cross_section[iE] = 0
        end
    end

    return cross_section ./ 1e4
end

function e_Oionion(Ep)
    # e_Oionion - O electron double ionization cross section (m^2)
    # Ep - electron energy (eV)
    # Source of data: Tima Sergienko, private communication
    cross_section = similar(Ep)

    for iE in eachindex(Ep)
        if Ep[iE] > 28.5 && Ep[iE] <= 250
            cross_section[iE] = (1-28.5/Ep[iE])*exp(-8.790669-22.50029*log(Ep[iE])+6.626340*(log(Ep[iE]))^2-0.8663901*(log(Ep[iE]))^3+0.04147531*(log(Ep[iE]))^4)
        elseif Ep[iE] > 250
            cross_section[iE] = 2.471376e-15*log(0.032*Ep[iE])/Ep[iE]
        else
            cross_section[iE] = 0
        end
    end

    return cross_section ./ 1e4
end
