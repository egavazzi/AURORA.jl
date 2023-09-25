using PyCall

function e_N2elastic(E) ✅
    cross_section = Vector{Float64}(undef, length(E))
    for ie in eachindex(E)
        if (E[ie] >= 1) & (E[ie] .< 1.733)
            cross_section[ie] = exp(-34.56132 + 0.9475102 * log(E[ie]) - 4.235151 * log(E[ie])^2 + 6.523666 * log(E[ie])^3)
        elseif (E[ie] >= 1.733) & (E[ie] .< 2.006)
            cross_section[ie] = exp(9.113062 - 242.3374 * log(E[ie]) + 441.5981 * log(E[ie])^2 - 262.3439 * log(E[ie])^3)
        elseif (E[ie] >= 2.006) & (E[ie] .< 2.317)
            cross_section[ie] = exp(72.26128 - 455.1942 * log(E[ie]) + 642.8157 * log(E[ie])^2 - 299.3418 * log(E[ie])^3)
        elseif (E[ie] >= 2.317) & (E[ie] .< 2.554)
            cross_section[ie] = exp(1767.169 - 6203.92 * log(E[ie]) + 7114.73 * log(E[ie])^2 - 2716.333 * log(E[ie])^3)
        elseif (E[ie] >= 2.554) & (E[ie] .< 2.775)
            cross_section[ie] = exp(867.3897 - 2850.537 * log(E[ie]) + 3000.209 * log(E[ie])^2 - 1050.929 * log(E[ie])^3)
        elseif (E[ie] >= 2.775) & (E[ie] .< 2.981)
            cross_section[ie] = exp(679.3589 - 2066.956 * log(E[ie]) + 1995.594 * log(E[ie])^2 - 641.9862 * log(E[ie])^3)
        elseif (E[ie] >= 2.981) & (E[ie] .< 3.215)
            cross_section[ie] = exp(-400.2765 + 865.4071 * log(E[ie]) - 668.1613 * log(E[ie])^2 + 167.3741 * log(E[ie])^3)
        elseif (E[ie] >= 3.215) & (E[ie] .< 3.457)
            cross_section[ie] = exp(1815.893 - 4716.509 * log(E[ie]) + 4004.182 * log(E[ie])^2 - 1132.109 * log(E[ie])^3)
        elseif (E[ie] >= 3.457) & (E[ie] .< 4.33)
            cross_section[ie] = exp(45.31069 - 184.5681 * log(E[ie]) + 142.8378 * log(E[ie])^2 - 36.88586 * log(E[ie])^3)
        elseif (E[ie] >= 4.33) & (E[ie] .< 1000)
            cross_section[ie] = exp(-35.26294 + 0.8019902 * log(E[ie]) - 0.202546 * log(E[ie])^2 + 0.007484 * log(E[ie])^3)
        elseif E[ie] >= 1000
            cross_section[ie] = 9.234e-14 / E[ie]
        else
            cross_section[ie] = 0
        end
    end
    cross_section .= cross_section ./ 1e4
    return cross_section
end

function e_N2rot0_2(E) ✅
    log10E = [-1.529205842868462, -1.401759744308854, -1.298455005566862,
        -1.229836683674959, -1.107192298263224, -0.930884014999626,
        -0.739165582563789, -0.534726081040409, -0.361197419278502,
        -0.210753800802626, -0.023659016727907, 0.148387783842545,
        0.195027357403972, 0.256843155275395, 0.271309779815279,
        0.284961440170103, 0.292280091257079, 0.313103310290649,
        0.339325032538764, 0.360990712227029, 0.383915274219258,
        0.406241931082478, 0.426693754073074, 0.438007445340180,
        0.442409564564855, 0.449441512598424, 0.461625748593174,
        0.475071834892396, 0.479153568268782, 0.698970004336019,
        1.000000000000000, 1.301029995663981, 1.477121254719663,
        1.698970004336019]
    log10Xs = [-16.489672666179871, -16.255790604873923, -16.087003295315228,
        -15.996264950740080, -15.953203452818123, -15.989191402537488,
        -16.041042208514060, -16.116569595034353, -16.168214419287263,
        -16.193450575738691, -16.195569244896330, -16.150452456808626,
        -16.107709885163285, -16.025430820560523, -15.961620468316118,
        -15.840488804947917, -15.803011086935593, -15.948567729023981,
        -15.533186544559602, -15.883157959631777, -15.476897837253116,
        -15.795237176295828, -15.528433895057361, -15.719342640373487,
        -15.765657146956446, -15.706366402956451, -15.622036628514593,
        -15.639548788652569, -15.655341886713167, -15.511449283499555,
        -15.414539270491499, -15.468521082957746, -15.552841968657781,
        -15.787812395596042]

    pyinterpolate = pyimport_conda("scipy.interpolate", "scipy");
    cross_section = 10 .^ pyinterpolate.PchipInterpolator(log10E, log10Xs)(log10.(E));
    cross_section = 10 .^ [pyinterpolate.PchipInterpolator(log10E, log10Xs)(log10.(E[E .<= 10^log10E[end-1]]))
                pyinterpolate.interp1d(log10E[end-1:end], log10Xs[end-1:end], kind="linear", fill_value="extrapolate")(log10.(E[E .>= 10^log10E[end-1]]))]
    I = findall(.!isfinite.(cross_section))
    cross_section[I] .= 0

    cross_section = cross_section / 1e4

    cross_section[E .> 10] .= 0 #TODO: FIX THIS/BG20190312  # wonder why? /EG20230924
    cross_section[E .< 0.001480105560000] .= 0

    return cross_section
end
