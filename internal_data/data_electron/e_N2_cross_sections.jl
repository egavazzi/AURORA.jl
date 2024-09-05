using PythonCall

"TO DO
elastic     ✅
rot0_2      ✅
rot0_4      ✅
rot0_6      ✅
rot0_8      ✅
vib0_1      ✅
vib0_2      ✅
vib0_3      ✅
vib0_4      ✅
vib0_5      ✅
vib0_6      ✅
vib0_7      ✅
a3sup       ✅
b3pg        ✅
w3du        ✅
bp3sum      ✅
ap1sum      ✅
w1du        ✅
e3sgp       ✅
ab1sgp      ✅
a1pg        ✅
c3pu        ✅
bp1sup      ✅
cp1sup      ✅
cp3pu       ✅
d3sup       ✅
f3pu        ✅
g3pu        ✅
M1M2        ✅
o1pu        ✅
dissociation      ✅
ionx2sgp    ✅
iona2pu     ✅
ionb2sup    ✅
dionv       ✅
ddion       ✅
"

function e_N2elastic(Ep)
    cross_section = Vector{Float64}(undef, length(Ep))
    for ie in eachindex(Ep)
        if (Ep[ie] >= 1) & (Ep[ie] .< 1.733)
            cross_section[ie] = exp(-34.56132 + 0.9475102 * log(Ep[ie]) - 4.235151 * log(Ep[ie])^2 + 6.523666 * log(Ep[ie])^3)
        elseif (Ep[ie] >= 1.733) & (Ep[ie] .< 2.006)
            cross_section[ie] = exp(9.113062 - 242.3374 * log(Ep[ie]) + 441.5981 * log(Ep[ie])^2 - 262.3439 * log(Ep[ie])^3)
        elseif (Ep[ie] >= 2.006) & (Ep[ie] .< 2.317)
            cross_section[ie] = exp(72.26128 - 455.1942 * log(Ep[ie]) + 642.8157 * log(Ep[ie])^2 - 299.3418 * log(Ep[ie])^3)
        elseif (Ep[ie] >= 2.317) & (Ep[ie] .< 2.554)
            cross_section[ie] = exp(1767.169 - 6203.92 * log(Ep[ie]) + 7114.73 * log(Ep[ie])^2 - 2716.333 * log(Ep[ie])^3)
        elseif (Ep[ie] >= 2.554) & (Ep[ie] .< 2.775)
            cross_section[ie] = exp(867.3897 - 2850.537 * log(Ep[ie]) + 3000.209 * log(Ep[ie])^2 - 1050.929 * log(Ep[ie])^3)
        elseif (Ep[ie] >= 2.775) & (Ep[ie] .< 2.981)
            cross_section[ie] = exp(679.3589 - 2066.956 * log(Ep[ie]) + 1995.594 * log(Ep[ie])^2 - 641.9862 * log(Ep[ie])^3)
        elseif (Ep[ie] >= 2.981) & (Ep[ie] .< 3.215)
            cross_section[ie] = exp(-400.2765 + 865.4071 * log(Ep[ie]) - 668.1613 * log(Ep[ie])^2 + 167.3741 * log(Ep[ie])^3)
        elseif (Ep[ie] >= 3.215) & (Ep[ie] .< 3.457)
            cross_section[ie] = exp(1815.893 - 4716.509 * log(Ep[ie]) + 4004.182 * log(Ep[ie])^2 - 1132.109 * log(Ep[ie])^3)
        elseif (Ep[ie] >= 3.457) & (Ep[ie] .< 4.33)
            cross_section[ie] = exp(45.31069 - 184.5681 * log(Ep[ie]) + 142.8378 * log(Ep[ie])^2 - 36.88586 * log(Ep[ie])^3)
        elseif (Ep[ie] >= 4.33) & (Ep[ie] .< 1000)
            cross_section[ie] = exp(-35.26294 + 0.8019902 * log(Ep[ie]) - 0.202546 * log(Ep[ie])^2 + 0.007484 * log(Ep[ie])^3)
        elseif Ep[ie] >= 1000
            cross_section[ie] = 9.234e-14 / Ep[ie]
        else
            cross_section[ie] = 0
        end
    end
    cross_section .= cross_section ./ 1e4
    return cross_section
end

function e_N2rot0_2(Ep)
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
    # import interpolate function from python
    pyinterpolate = pyimport("scipy.interpolate")
    # cross_section = 10 .^ pyinterpolate.PchipInterpolator(log10E, log10Xs)(log10.(Ep));
    cross_section = 10 .^ [pyinterpolate.PchipInterpolator(log10E, log10Xs)(log10.(Ep[Ep .<= 10^log10E[end-1]])),
                pyinterpolate.interp1d(log10E[end-1:end], log10Xs[end-1:end], kind="linear", fill_value="extrapolate")(log10.(Ep[Ep .>= 10^log10E[end-1]]))]
    # convert from a Python array back to a Julia array
    cross_section = pyconvert.(Array, cross_section)
    # merge the two parts
    cross_section = vcat(cross_section[1], cross_section[2])

    I = findall(.!isfinite.(cross_section))
    cross_section[I] .= 0

    cross_section = cross_section / 1e4

    cross_section[Ep .> 10] .= 0 #TODO: FIX THIS/BG20190312  # wonder why? /EG20230924
    cross_section[Ep .< 0.001480105560000] .= 0

    return cross_section
end

function e_N2rot0_4(Ep)
    log10E = [-1.532680902603447, -1.404874076776330, -1.275825566920618,
              -1.224367870665778, -1.149260425122097, -0.985947211351686,
              -0.834892017817066, -0.642148004931654, -0.267182529393861,
              -0.087071929590095, -0.002885752451284, 0.091967356784349,
              0.160395697309093, 0.195408146635714, 0.219793807440825,
              0.246011772608117, 0.282549189711953, 0.290486080963258,
              0.293629539250526, 0.298637983878642, 0.305483198238818,
              0.311674551718697, 0.316057534684070, 0.322827300608398,
              0.327398396418622, 0.335001803552344, 0.340464631369559,
              0.347584582685744, 0.361414574811236, 0.368229005640144,
              0.377677383680488, 0.380158540333515, 0.382076382389718,
              0.384787439246178, 0.397755986770893, 0.401409239603281,
              0.406827163211131, 0.409948376785994, 0.418183570165017,
              0.422882192428233, 0.426271155084866, 0.428770588446716,
              0.440721896378690, 0.443789449978017, 0.447612136905164,
              0.460363315178215, 0.464043471225615, 0.479034714667705,
              0.698970004336019, 1.000000000000000, 1.301029995663981,
              1.477121254719663, 1.698970004336019]

    log10Xs = [-17.306560217472983, -16.905296630034872, -16.496472163619632,
               -16.346950305478465, -16.302769455158653, -16.372165876080111,
               -16.426413029943376, -16.531167125488427, -16.790591014528083,
               -16.912714248287095, -16.943688489526014, -16.912239953864546,
               -16.800436413819938, -16.590216924422101, -16.405144814674895,
               -16.102766217156070, -15.589673466968044, -15.508545593737043,
               -15.500300448742312, -15.580114743299788, -15.801601894983513,
               -15.947592289811023, -15.842459615546449, -15.608113401394665,
               -15.373307955548082, -15.245568441986372, -15.205532305410159,
               -15.427135164359850, -15.943514591785419, -15.680492955101307,
               -15.223843323199201, -15.176424582771727, -15.160459552129517,
               -15.186373523082741, -15.698794908820723, -15.795806879885383,
               -15.851026564276077, -15.807316835052490, -15.432362387183495,
               -15.326658892601827, -15.300812117955953, -15.324945790826710,
               -15.784820589642830, -15.844461318346980, -15.785259390571994,
               -15.519190930094789, -15.490189382908351, -15.625025350388160,
               -15.804100347590765, -15.931814138253838, -16.080921907623924,
               -15.903089986991944, -15.860120913598763]

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = 10 .^ pyinterpolate.PchipInterpolator(log10E, log10Xs)(log10.(Ep));
    cross_section = pyconvert.(Array, cross_section)
    I = findall(.!isfinite.(cross_section))
    cross_section[I] .= 0

    cross_section = cross_section / 1e4

    cross_section[Ep .> 10] .= 0 #TODO: FIX THIS/BG20190312  # wonder why? /EG20230924
    cross_section[Ep .< 0.004933884000000] .= 0

    return cross_section
end

function e_N2rot0_6(Ep)
    log10E = [-1.540835374081948, -1.393951276348084, -1.299408565019692,
              -1.198041825677987, -1.085443481487013, -1.027862686921294,
              -0.853878202151661, -0.718632215072602, -0.500218151334885,
              -0.387657446570351, -0.358603573302421, 1.000000000000000,
              1.301029995663981, 1.477121254719663, 1.698970004336019]

    log10Xs = [-17.463274018565434, -16.993190004022189, -16.664905513524634,
               -16.587790797622958, -16.660806579985472, -16.726677458224287,
               -16.816693087537782, -16.893226274300982, -17.065268564666383,
               -17.141544862344141, -17.151656694256861, -18.154901959985743,
               -17.096910013008056, -16.443697499232712, -16.585026652029178]

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = 10 .^ pyinterpolate.PchipInterpolator(log10E, log10Xs)(log10.(Ep));
    cross_section = pyconvert.(Array, cross_section)
    I = findall(.!isfinite.(cross_section))
    cross_section[I] .= 0

    cross_section = cross_section / 1e4

    cross_section[Ep .> 10] .= 0 #TODO: FIX THIS/BG20190312  # wonder why? /EG20230924
    cross_section[Ep .< 0.010361812440000] .= 0

    return cross_section
end

function e_N2rot0_8(Ep)
    log10E = [-1.528601289020590, -1.426292957506456, -1.329246907900923,
              -1.276205819361377, -1.202716630014793, -1.121888884265058,
              -1.049473630001863, -0.957742617010783, -0.856184155958383,
              -0.785278063101844, -0.658159947064807, -0.576218700202296,
              -0.413581594354676]

    log10Xs = [-17.416103498733147, -17.126710494427051, -16.871536453251654,
               -16.774198368667303, -16.732590449548407, -16.757990427343014,
               -16.809412134822608, -16.874111071306828, -16.930759733883072,
               -16.972370479146925, -17.063587133755963, -17.133072146723649,
               -17.275381828471725]

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = 10 .^ pyinterpolate.PchipInterpolator(log10E, log10Xs)(log10.(Ep));
    cross_section = pyconvert.(Array, cross_section)
    I = findall(.!isfinite.(cross_section))
    cross_section[I] .= 0

    cross_section = cross_section / 1e4

    cross_section[Ep .> 10] .= 0 #TODO: FIX THIS/BG20190312  # wonder why? /EG20230924
    cross_section[Ep .< 0.017764640640000] .= 0

    return cross_section
end

function e_N2vib0_1(Ep)
    E = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    s = [1e-23, 3e-23, 3.8e-23, 4.2e-23, 4.6e-23, 5.5e-23, 6.7e-23, 8.1e-23]

    E = vcat(E, [1.6, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2, 2.05, 2.1, 2.15, 2.2, 2.23, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.73, 2.8, 2.85, 2.9, 3, 3.05, 3.075, 3.14, 3.2, 3.3, 3.4])
    s = vcat(s, [0.28e-20, 0.51e-20, 1.2e-20, 1.9e-20, 2.6e-20, 4.5e-20, 5.7e-20, 4.2e-20, 2.2e-20, 1.4e-20, 2.5e-20, 5.0e-20, 5.8e-20, 4.4e-20, 3.6e-20, 1.9e-20, 1.5e-20, 2.5e-20, 4.5e-20, 3.6e-20, 1.67e-20, 1.21e-20, 2.78e-20, 3.20e-20, 2.04e-20, 0.88e-20, 1.35e-20, 1.92e-20, 1.58e-20, 0.79e-20, 0.97e-20, 0.56e-20])

    E = vcat(E, [5, 7.5, 10, 15, 18, 20, 23, 25, 30, 50, 75])
    s = vcat(s, [6.5e-22, 3.1e-22, 1.4e-22, 4.1e-22, 7.5e-22, 19.4e-22, 12.1e-22, 7.2e-22, 2.4e-22, 1.4e-22, 0.67e-22])

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = [pyinterpolate.PchipInterpolator(E, log.(s))(Ep[Ep .< E[end]]),
               pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    cross_section = pyconvert.(Array, cross_section)
    cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = exp.(cross_section)

    I = findall(.!isfinite.(cross_section))
    cross_section[I] .= 0

    cross_section[Ep .> 10] .= 0 #TODO: FIX THIS/BG20190312  # wonder why? /EG20230924
    cross_section[Ep .< 0.2888] .= 0

    return cross_section
end

function e_N2vib0_2(Ep)
    E = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    s = [1, 3, 3.8, 4.2, 4.6, 5.5, 6.7, 8.1] .* 1e-23
    E = vcat(E, [1.67, 1.75, 1.8, 1.85, 1.9, 1.9375, 1.978, 2.0, 2.05, 2.1, 2.125, 2.2, 2.24, 2.28, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4])
    s = vcat(s, [0.08, 0.19367, 0.28, 1.37, 1.94, 3.82, 4.04, 3.65, 3.16, 2.24, 1.79, 0.57, 1.14, 2.47, 4.16, 3.88, 2.16, 1.25, 0.51, 0.97, 1.94, 2.32, 1.67, 0.97, 0.4, 0.51, 0.91, 1.2, 1.08, 0.39, 0.23, 0.54, 0.74, 0.46, 0.17, 0.20] .* 1e-20)
    E = vcat(E, [5, 7.5, 10, 15, 18, 20, 23, 25, 30, 50, 75])
    s = vcat(s, [6.5, 3.1, 1.4, 4.1, 7.5, 19.4, 12.1, 7.2, 2.4, 1.4, 0.67] .* 1e-22)

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = [pyinterpolate.PchipInterpolator(E, log.(s))(Ep[Ep .< E[end]]),
                     pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    cross_section = pyconvert.(Array, cross_section)
    cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = exp.(cross_section)

    I = findall(.!isfinite.(cross_section))
    cross_section[I] .= 0

    cross_section[Ep .> 10] .= 0 #TODO: FIX THIS/BG20190312  # wonder why? /EG20240515
    cross_section[Ep .< 0.2888] .= 0

    return cross_section
end

function e_N2vib0_3(Ep)
    E = [1.0271, 1.8119, 1.9281, 2.3916, 2.4163, 2.5739, 2.7995, 3.0116, 3.1689, 3.3487, 3.5064, 3.7333, 3.8524, 4.0669, 4.067, 4.2041, 4.3523, 4.5029, 4.6469]
    s = [7.0427e-25, 5.5034e-24, 2.5652e-22, 1.673e-20, 1.6872e-20, 3.4988e-21, 1.5817e-20, 2.1678e-21, 1.1187e-20, 1.5309e-21, 4.3424e-21, 1.2751e-21, 1.7867e-21, 7.6565e-22, 7.6526e-22, 1.0711e-21, 5.0849e-22, 1.0811e-21, 5.5485e-22]

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = pyinterpolate.PchipInterpolator(E, s)(Ep)
    cross_section = pyconvert.(Array, cross_section)
    cross_section = cross_section .* (Ep .< E[end])
    # Handling energies beyond the last data point using e_N2vib0_1
    cross_section_01 = e_N2vib0_1(Ep) .* (Ep .> E[end])
    cross_section = cross_section .+ cross_section_01

    cross_section[Ep .> 10] .= 0 #TODO: FIX THIS/BG20190312
    cross_section[Ep .< 0.8559] .= 0

    return cross_section
end

function e_N2vib0_4(Ep)
    # Initial combination of e_N2vib0_1 and e_N2vib0_2
    # initial_cross_section = (0.3 * e_N2vib0_1(Ep) + 0.35 * e_N2vib0_2(Ep)) * 2 / 5 # not used? /EG20240515

    # Energy and cross-section data
    E = [1.6667, 1.8867, 2.2155, 2.535, 2.7699, 3.049, 3.1412, 3.2972, 3.3004, 3.4765,  3.4881, 3.6861, 3.8694, 4.0781]
    s = [1.9526e-23, 5.5125e-23, 1.1669e-20, 7.5773e-22, 1.0646e-20, 4.0738e-21, 4.0361e-21, 2.314e-21, 2.293e-21, 2.2845e-21, 2.2741e-21, 1.2528e-21, 9.8467e-22, 5.0127e-22]
    E = vcat([1.1342], E)
    s = vcat([1e-27], s)

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = pyinterpolate.PchipInterpolator(E, s)(Ep)
    cross_section = pyconvert.(Array, cross_section)
    cross_section = cross_section .* (Ep .< E[end])

    # Handling energies beyond the last data point using e_N2vib0_1
    cross_section_01 = e_N2vib0_1(Ep) .* (Ep .> E[end])
    cross_section = cross_section .+ 0.35378 .* cross_section_01

    # Set cross-section to zero outside specified energy ranges
    cross_section[Ep .> 10] .= 0 #TODO: FIX THIS/BG20190312
    cross_section[Ep .< 1.1342] .= 0

    return cross_section
end

function e_N2vib0_5(Ep)
    E = [1.7459, 2.0816, 2.3099, 2.4942, 2.6043, 2.8488, 3.1042, 3.2912, 3.4976, 3.7519, 3.8671]
    s = [6.8315e-24, 2.1617e-23, 2.9028e-21, 3.9878e-21, 6.2584e-21, 2.6986e-22, 3.9496e-21, 2.6388e-22, 1.613e-21, 1.7546e-22, 5.1621e-22]
    s = vcat([s[1]/3], s)
    E = vcat([1.4088], E)

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = pyinterpolate.PchipInterpolator(E, s)(Ep)
    cross_section = pyconvert.(Array, cross_section)
    cross_section = cross_section .* (Ep .< E[end])

    cross_section_01 = e_N2vib0_1(Ep) .* (Ep .> E[end])
    cross_section = cross_section .+ cross_section_01 / 10.663 * 2.5

    cross_section[Ep .> 10] .= 0  # TODO: FIX THIS/BG20190312
    cross_section[Ep .< 1.4088] .= 0

    return cross_section
end

function e_N2vib0_6(Ep)
    E = [1.5176, 2.143, 2.2964, 2.5275, 2.8863, 3.0358, 3.1708, 3.4811]
    s = [1.0846e-25, 9.2113e-24, 4.9432e-23, 4.459e-21, 3.4105e-22, 1.1193e-21, 2.9207e-22, 2.1537e-22]

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = pyinterpolate.PchipInterpolator(E, s)(Ep)
    cross_section = pyconvert.(Array, cross_section)
    cross_section = cross_section .* (Ep .< E[end])

    cross_section_01 = e_N2vib0_1(Ep) .* (Ep .> E[end])
    cross_section = cross_section .+ cross_section_01 / 20

    cross_section[Ep .> 10] .= 0  # TODO: FIX THIS/BG20190312
    cross_section[Ep .< 1.68] .= 0

    return cross_section
end

function e_N2vib0_7(Ep)
    E = [2.0073, 2.2342, 2.3663, 2.5946, 2.8222, 3.0483, 3.2516, 3.3851, 3.5291, 4.1979]
    s = [6.317e-23, 1.3348e-22, 6.3896e-22, 2.4423e-21, 4.2059e-22, 1.6918e-21, 3.8452e-22, 7.5752e-22, 3.9773e-22, 2.5062e-22]

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = pyinterpolate.PchipInterpolator(E, s)(Ep)
    cross_section = pyconvert.(Array, cross_section)
    cross_section = cross_section .* (Ep .< E[end])

    cross_section_01 = e_N2vib0_1(Ep) .* (Ep .> E[end]) / 5
    cross_section = cross_section .+ cross_section_01

    cross_section[Ep .> 10] .= 0  # TODO: FIX THIS/BG20190312
    cross_section[Ep .< 1.9] .= 0

    return cross_section
end

function e_N2a3sup(Ep)
    E = [7, 8, 9, 10, 13, 15, 17, 20, 30, 50]
    s = [.2, 2, 4.5, 15, 17.6, 20.5, 22, 15, 6, 3.8] .* 1e-22

    Xs = zeros(length(Ep))

    for ie in eachindex(Ep)
        if Ep[ie] > 6.17 && Ep[ie] < 11.096
            Xs[ie] = (1 - 6.17 / Ep[ie]) * exp(-48.28302 - 501.0230 * log(Ep[ie]) + 635.6719 * log(Ep[ie])^2 - 264.7979 * log(Ep[ie])^3 + 36.53586 * log(Ep[ie])^4)
        elseif Ep[ie] >= 11.096 && Ep[ie] < 17.05
            Xs[ie] = exp(4550.309 - 6733.425 * log(Ep[ie]) + 3697.708 * log(Ep[ie])^2 - 900.7303 * log(Ep[ie])^3 + 82.11416 * log(Ep[ie])^4)
        elseif Ep[ie] >= 17.05 && Ep[ie] < 100.0
            Xs[ie] = exp(-142.1066 + 124.6773 * log(Ep[ie]) - 55.47813 * log(Ep[ie])^2 + 10.80510 * log(Ep[ie])^3 - 0.7832104 * log(Ep[ie])^4)
        elseif Ep[ie] >= 100
            Xs[ie] = 9.608389e-13 / Ep[ie]^3
        else
            Xs[ie] = 0
        end
    end

    pyinterpolate = pyimport("scipy.interpolate")
    XsItikawa = [pyinterpolate.PchipInterpolator(E, log.(s))(Ep[Ep .< E[end]]),
                 pyinterpolate.interp1d(log.(E), log.(s), fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    XsItikawa = pyconvert.(Array, XsItikawa)
    XsItikawa = vcat(XsItikawa[1], XsItikawa[2])
    XsItikawa = exp.(XsItikawa)

    XsItikawa[.!isfinite.(XsItikawa)] .= 0

    Xs = [XsItikawa[Ep .<= 25]; Xs[Ep .> 25] / 1e4]
    Xs[Ep .< 6.1688] .= 0

    return Xs
end

function e_N2b3pg(Ep)
    E = [7, 8, 9, 10, 12.5, 15, 17, 20, 30, 50]
    s = [.2, 2, 4.5, 25, 32, 22.5, 18.5, 13, 8, 3] .* 1e-22

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = [pyinterpolate.PchipInterpolator(E, log.(s))(Ep[Ep .< E[end]]),
                     pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    cross_section = pyconvert.(Array, cross_section)
    cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section_temp = similar(Ep)
    for ie in eachindex(Ep)
        if Ep[ie] > 7.35 && Ep[ie] < 15.52
            cross_section_temp[ie] = (1 - 7.35 / Ep[ie]) * exp(-1139.542 + 1583.892 * log(Ep[ie]) - 844.6222 * log(Ep[ie])^2 + 198.2095 * log(Ep[ie])^3 - 17.29356 * log(Ep[ie])^4)
        elseif Ep[ie] >= 15.52 && Ep[ie] < 43.1
            cross_section_temp[ie] = exp(-0.9660707 - 18.27795 * log(Ep[ie]) - 5.686039 * log(Ep[ie])^2 + 4.259477 * log(Ep[ie])^3 - 0.5800993 * log(Ep[ie])^4)
        elseif Ep[ie] >= 43.1 && Ep[ie] < 100.0
            cross_section_temp[ie] = exp(-1087.040 + 1098.993 * log(Ep[ie]) - 429.7112 * log(Ep[ie])^2 + 74.32909 * log(Ep[ie])^3 - 4.807258 * log(Ep[ie])^4)
        elseif Ep[ie] >= 100
            cross_section_temp[ie] = 6.148782e-13 / Ep[ie]^3
        else
            cross_section_temp[ie] = 0
        end
    end

    cross_section[Ep .> 95] = cross_section_temp[Ep .> 95] / 1e4
    cross_section[Ep .< 7.3532] .= 0

    return cross_section
end

function e_N2w3du(Ep)
    cross_section_temp = similar(Ep)

    for ie in eachindex(Ep)
        if Ep[ie] > 7.36 && Ep[ie] < 100
            cross_section_temp[ie] = (1 - 7.36 / Ep[ie]) * exp(-275.1462 + 274.0435 * log(Ep[ie]) - 116.3761 * log(Ep[ie])^2 + 21.57467 * log(Ep[ie])^3 - 1.484377 * log(Ep[ie])^4)
        elseif Ep[ie] >= 100
            cross_section_temp[ie] = 7.036322e-13 / Ep[ie]^3
        else
            cross_section_temp[ie] = 0
        end
    end

    E = [7, 8, 9, 10, 12.5, 15, 17, 20, 30, 50]
    s = [.2, 2, 4.5, 7.5, 24.5, 33, 34, 31, 7, 2.2] .* 1e-22

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = [pyinterpolate.PchipInterpolator(E, log.(s))(Ep[Ep .< E[end]]),
                     pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    cross_section = pyconvert.(Array, cross_section)
    cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 7.3622] .= 0
    cross_section[Ep .> 28.25] = cross_section_temp[Ep .> 28.25] / 1e4

    return cross_section
end

function e_N2bp3sum(Ep)
    cross_section_temp = similar(Ep)

    for ie in eachindex(Ep)
        if Ep[ie] > 8.16 && Ep[ie] < 29.55
            cross_section_temp[ie] = (1 - 8.16 / Ep[ie]) * exp(-667.3893 + 764.3953 * log(Ep[ie]) - 346.0872 * log(Ep[ie])^2 + 69.19737 * log(Ep[ie])^3 - 5.167390 * log(Ep[ie])^4)
        elseif Ep[ie] >= 29.55 && Ep[ie] < 100
            cross_section_temp[ie] = exp(-34.27635 + 0.0 * log(Ep[ie]) - 1.745869 * log(Ep[ie])^2 + 0.5473512 * log(Ep[ie])^3 - 0.0553294 * log(Ep[ie])^4)
        elseif Ep[ie] >= 100
            cross_section_temp[ie] = 2.77e-13 / Ep[ie]^3
        else
            cross_section_temp[ie] = 0
        end
    end

    E = [8.1647, 12, 15, 18, 20, 30, 50]
    s = [.0001, 8, 13, 8, 4, 3, 2] .* 1e-22

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = [pyinterpolate.PchipInterpolator(E, log.(s))(Ep[Ep .< E[end]]),
                     pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    cross_section = pyconvert.(Array, cross_section)
    cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 8.1647] .= 0
    cross_section[Ep .> 30] = cross_section_temp[Ep .> 30] / 1e4

    return cross_section
end

function e_N2ap1sum(Ep)
    cross_section_temp = similar(Ep)

    for ie in eachindex(Ep)
        if Ep[ie] > 8.4 && Ep[ie] < 29.937
            cross_section_temp[ie] = (1 - 8.4 / Ep[ie]) * exp(-2652.316 + 3473.205 * log(Ep[ie]) - 1723.654 * log(Ep[ie])^2 + 378.6041 * log(Ep[ie])^3 - 31.07141 * log(Ep[ie])^4)
        elseif Ep[ie] >= 29.937 && Ep[ie] < 100
            cross_section_temp[ie] = exp(-436.8312 + 412.7602 * log(Ep[ie]) - 159.5703 * log(Ep[ie])^2 + 27.18180 * log(Ep[ie])^3 - 1.728542 * log(Ep[ie])^4)
        elseif Ep[ie] >= 100
            cross_section_temp[ie] = 2.495837e-17 / Ep[ie]
        else
            cross_section_temp[ie] = 0
        end
    end

    E = [8.3987, 15, 18, 20, 30, 50]
    s = [.0001, 11, 5, 3.5, 2, .8] .* 1e-22

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = [pyinterpolate.PchipInterpolator(E, log.(s))(Ep[Ep .< E[end]]),
                     pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    cross_section = pyconvert.(Array, cross_section)
    cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 8.3987] .= 0
    cross_section[Ep .> 45] = cross_section_temp[Ep .> 45] / 1e4

    return cross_section
end

function e_N2w1du(Ep)
    cross_section_temp = similar(Ep)

    for ie in eachindex(Ep)
        if Ep[ie] > 8.89 && Ep[ie] < 20.487
            cross_section_temp[ie] = (1 - 8.89 / Ep[ie]) * exp(-5231.492 + 7031.762 * log(Ep[ie]) - 3566.315 * log(Ep[ie])^2 + 802.9877 * log(Ep[ie])^3 - 67.74330 * log(Ep[ie])^4)
        elseif Ep[ie] >= 20.487 && Ep[ie] < 100
            cross_section_temp[ie] = exp(-131.5858 + 104.4165 * log(Ep[ie]) - 43.27655 * log(Ep[ie])^2 + 7.735902 * log(Ep[ie])^3 - 0.5085983 * log(Ep[ie])^4)
        elseif Ep[ie] >= 100
            cross_section_temp[ie] = 7.415168e-17 / Ep[ie]
        else
            cross_section_temp[ie] = 0
        end
    end

    E = [8.8895, 12, 15, 17, 20, 30, 50]
    s = [.0001, 12, 9.5, 7, 3, 1.75, 0.7] .* 1e-22

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = [pyinterpolate.PchipInterpolator(E, log.(s))(Ep[Ep .< E[end]]),
                     pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    cross_section = pyconvert.(Array, cross_section)
    cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 8.8895] .= 0
    cross_section[30 .<= Ep .< 35] = (cross_section[30 .<= Ep .< 35] * 2 / 3 + cross_section_temp[30 .<= Ep .< 35] * 1 / 3 / 1e4)
    cross_section[35 .<= Ep .< 40] = (cross_section[35 .<= Ep .< 40] * 2 / 3 + cross_section_temp[35 .<= Ep .< 40] * 1 / 3 / 1e4)
    cross_section[40 .<= Ep] = cross_section_temp[40 .<= Ep] / 1e4

    return cross_section
end

function e_N2e3sgp(Ep)
    cross_section = zeros(length(Ep))

    for ie in eachindex(Ep)
        if Ep[ie] > 11.88 && Ep[ie] < 30.35
            cross_section[ie] = (1 - 11.88 / Ep[ie]) * exp(74.2133 - 124.664 * log(Ep[ie]) + 44.6708 * log(Ep[ie])^2 - 5.30151 * log(Ep[ie])^3)
        elseif Ep[ie] >= 30.35 && Ep[ie] < 50
            cross_section[ie] = exp(-105.2329 + 54.4886 * log(Ep[ie]) - 14.77986 * log(Ep[ie])^2 + 1.23869 * log(Ep[ie])^3)
        elseif Ep[ie] >= 50
            cross_section[ie] = 8.7619e-15 / Ep[ie]^3
        else
            cross_section[ie] = 0
        end
    end

    cross_section = cross_section / 1e4

    return cross_section
end

function e_N2ab1sgp(Ep)
    cross_section = zeros(length(Ep))

    for ie in eachindex(Ep)
        if Ep[ie] > 12.25 && Ep[ie] < 24.98
            cross_section[ie] = (1 - 12.25 / Ep[ie]) * exp(91.77961 - 148.1616 * log(Ep[ie]) + 55.75255 * log(Ep[ie])^2 - 6.95604 * log(Ep[ie])^3)
        elseif Ep[ie] >= 24.98 && Ep[ie] < 50
            cross_section[ie] = exp(80.2784 - 92.0627 * log(Ep[ie]) + 23.36969 * log(Ep[ie])^2 - 1.985404 * log(Ep[ie])^3)
        elseif Ep[ie] >= 50
            cross_section[ie] = 7.143e-17 / Ep[ie]
        else
            cross_section[ie] = 0
        end
    end

    cross_section = cross_section / 1e4

    cross_section[.!isfinite.(cross_section)] .= 0

    return cross_section
end

function e_N2a1pg(Ep)
    cross_section_temp = zeros(length(Ep))

    for ie in eachindex(Ep)
        if Ep[ie] > 8.55 && Ep[ie] < 100.0
            cross_section_temp[ie] = (1 - 8.55 / Ep[ie]) * exp(-108.6546 + 73.07788 * log(Ep[ie]) - 27.46333 * log(Ep[ie])^2 + 4.465812 * log(Ep[ie])^3 - 0.2689957 * log(Ep[ie])^4)
        elseif Ep[ie] >= 100.0
            cross_section_temp[ie] = 7.207318e-16 / Ep[ie]
        else
            cross_section_temp[ie] = 0
        end
    end

    E = [10, 12.5, 15, 16.25, 17.5, 20, 22.0, 30.0, 35.0, 40.0, 50, 60, 75, 90, 100]
    s = [0.42, 2.29, 3.67, 4, 3.67, 3.08, 2.67, 2.12, 1.67, 1.5, 1.17, 0.83, 0.75, 0.57, 0.5] .* 1e-21

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = [pyinterpolate.PchipInterpolator(E, log.(s))(Ep[Ep .< E[end]]),
                     pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    cross_section = pyconvert.(Array, cross_section)
    cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 8.5489] .= 0
    cross_section[Ep .> 31.25] = cross_section_temp[Ep .> 31.25] / 1e4

    return cross_section
end

function e_N2c3pu(Ep)
    cross_section_temp = similar(Ep)

    for ie in eachindex(Ep)
        if Ep[ie] > 11.03 && Ep[ie] < 20.169
           cross_section_temp[ie] = (1 - 11.03 / Ep[ie]) * exp(-9134.460 + 12303.41 * log(Ep[ie]) - 6226.840 * log(Ep[ie])^2 + 1397.931 * log(Ep[ie])^3 - 117.4893 * log(Ep[ie])^4)
        elseif Ep[ie] >= 20.169 && Ep[ie] < 100
           cross_section_temp[ie] = exp(-145.6415 + 117.1985 * log(Ep[ie]) - 47.66838 * log(Ep[ie])^2 + 8.547379 * log(Ep[ie])^3 - 0.5779936 * log(Ep[ie])^4)
        elseif Ep[ie] >= 100
           cross_section_temp[ie] = 5.5409940e-13 / Ep[ie]^3
        else
           cross_section_temp[ie] = 0
        end
    end

    E = [12.1, 12.67, 13.33, 13.5, 15, 17, 20, 30, 50]
    s = [0.5, 1, 2, 3, 3.75, 2.05, 1.44, 0.54, 0.21] .* 1e-21

    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = [pyinterpolate.PchipInterpolator(E, log.(s))(Ep[Ep .< E[end]]),
                     pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    cross_section = pyconvert.(Array, cross_section)
    cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    cross_section[Ep .< 11.032] .= 0
    cross_section[Ep .> 24.4] = cross_section_temp[Ep .> 24.4] / 1e4

    return cross_section
end

function e_N2bp1sup(Ep)
    cross_section = similar(Ep)

    for iE in length(Ep):-1:1
        if Ep[iE] > 12.85 && Ep[iE] < 152
            cross_section[iE] = (1 - 12.85 / Ep[iE]) * exp(-42.893187 + 3.656131 * log(Ep[iE]) - 0.9243795 * log(Ep[iE])^2 + 0.0642402 * log(Ep[iE])^3)
        elseif Ep[iE] >= 152
            cross_section[iE] = 3.649e-16 * log(0.0573 * Ep[iE]) / Ep[iE]
        else
            cross_section[iE] = 0
        end
    end

    cross_section *= 1e-4

    return cross_section
end

function e_N2cp1sup(Ep)
    cross_section = similar(Ep)

    for iE in length(Ep):-1:1
        if Ep[iE] > 12.94 && Ep[iE] < 151
            cross_section[iE] = (1 - 12.94 / Ep[iE]) * exp(-42.69132 + 3.656131 * log(Ep[iE]) - 0.9243795 * log(Ep[iE])^2 + 0.0642402 * log(Ep[iE])^3)
        elseif Ep[iE] >= 151
            cross_section[iE] = 4.466e-16 * log(0.0573 * Ep[iE]) / Ep[iE]
        else
            cross_section[iE] = 0
        end
    end

    cross_section *= 1e-4

    return cross_section
end

function e_N2cp3pu(Ep)
    cross_section = similar(Ep)

    for iE in length(Ep):-1:1
        if Ep[iE] > 12.08 && Ep[iE] < 18.36
            cross_section[iE] = (1 - 12.08 / Ep[iE]) * exp(-3611.089 + 3701.604 * log(Ep[iE]) - 1276.938 * log(Ep[iE])^2 + 146.5513 * log(Ep[iE])^3)
        elseif Ep[iE] >= 18.36 && Ep[iE] < 90
            cross_section[iE] = exp(33.71613 - 61.7778 * log(Ep[iE]) + 16.6116 * log(Ep[iE])^2 - 1.50206 * log(Ep[iE])^3)
        elseif Ep[iE] >= 90
            cross_section[iE] = 2.616e-14 / Ep[iE]^3
        else
            cross_section[iE] = 0
        end
    end

    cross_section *= 1e-4

    return cross_section
end

function e_N2d3sup(Ep)
    cross_section = similar(Ep)

    for iE in length(Ep):-1:1
        if Ep[iE] > 12.85 && Ep[iE] < 50
            cross_section[iE] = (1 - 12.85 / Ep[iE]) * exp(0.938347 - 31.08899 * log(Ep[iE]) + 8.257418 * log(Ep[iE])^2 - 0.793736 * log(Ep[iE])^3)
        elseif Ep[iE] >= 50
            cross_section[iE] = 6.3166e-14 / Ep[iE]^3
        else
            cross_section[iE] = 0
        end
    end

    cross_section *= 1e-4

    return cross_section
end

function e_N2f3pu(Ep)
    cross_section = similar(Ep)

    for iE in length(Ep):-1:1
        if Ep[iE] > 12.75 && Ep[iE] < 23.28
            cross_section[iE] = (1 - 12.75 / Ep[iE]) * exp(27.4082 - 52.0658 * log(Ep[iE]) + 13.78086 * log(Ep[iE])^2 - 1.26079 * log(Ep[iE])^3)
        elseif Ep[iE] >= 23.28 && Ep[iE] < 90
            cross_section[iE] = exp(1.13172 - 32.235 * log(Ep[iE]) + 8.5621 * log(Ep[iE])^2 - 0.78723 * log(Ep[iE])^3)
        elseif Ep[iE] >= 90
            cross_section[iE] = 3.171e-13 / Ep[iE]^3
        else
            cross_section[iE] = 0
        end
    end

    cross_section *= 1e-4

    return cross_section
end

function e_N2g3pu(Ep)
    cross_section = similar(Ep)

    for iE in length(Ep):-1:1
        if Ep[iE] > 12.8 && Ep[iE] < 18.84
            cross_section[iE] = (1 - 12.8 / Ep[iE]) * exp(163.9467 - 189.6174 * log(Ep[iE]) + 60.6054 * log(Ep[iE])^2 - 6.615773 * log(Ep[iE])^3)
        elseif Ep[iE] >= 18.84 && Ep[iE] < 80
            cross_section[iE] = exp(6.881715 - 37.398 * log(Ep[iE]) + 10.292 * log(Ep[iE])^2 - 0.97713 * log(Ep[iE])^3)
        elseif Ep[iE] >= 80
            cross_section[iE] = 4.441e-13 / Ep[iE]^3
        else
            cross_section[iE] = 0
        end
    end

    cross_section *= 1e-4

    return cross_section
end

function e_N2M1M2(Ep)
    cross_section = similar(Ep)

    for iE in length(Ep):-1:1
        if Ep[iE] > 13.15 && Ep[iE] < 24.0
            cross_section[iE] = (1 - 13.15 / Ep[iE]) * exp(115.6489 - 142.1746 * log(Ep[iE]) + 44.0739 * log(Ep[iE])^2 - 4.637555 * log(Ep[iE])^3)
        elseif Ep[iE] >= 24.0 && Ep[iE] < 80
            cross_section[iE] = exp(5.57238 - 37.398 * log(Ep[iE]) + 10.292 * log(Ep[iE])^2 - 0.97713 * log(Ep[iE])^3)
        elseif Ep[iE] >= 80
            cross_section[iE] = 1.199e-13 / Ep[iE]^3
        else
            cross_section[iE] = 0
        end
    end

    cross_section *= 1e-4

    return cross_section
end

function e_N2o1pu(Ep)
    cross_section = similar(Ep)

    for iE in length(Ep):-1:1
        if Ep[iE] > 13.1 && Ep[iE] < 150
            cross_section[iE] = (1 - 13.1 / Ep[iE]) * exp(-44.40084 + 3.656131 * log(Ep[iE]) - 0.9243795 * log(Ep[iE])^2 + 0.0642402 * log(Ep[iE])^3)
        elseif Ep[iE] >= 150
            cross_section[iE] = 8.0809e-17 * log(0.0573 * Ep[iE]) / Ep[iE]
        else
            cross_section[iE] = 0
        end
    end

    cross_section *= 1e-4

    return cross_section
end

function e_N2dissociation(Ep)
    E = [20, 23.6, 27.8, 35, 45, 55, 65, 80, 96, 110, 125, 148, 170, 195, 245, 295, 2985.1, 49183, 55306]
    s = [0.87, 1.13, 1.39, 1.54, 1.7, 1.87, 1.87, 2.04, 2.07, 1.96, 1.9, 1.87, 1.78, 1.74, 1.57, 1.48, 5.0198e-1/3, 4.5953e-2/3, 4.47e-2/3] * 1e-20


    pyinterpolate = pyimport("scipy.interpolate")
    cross_section = [pyinterpolate.PchipInterpolator(E, log.(s))(Ep[Ep .< E[end]]),
                     pyinterpolate.interp1d(log.(E), log.(s), kind="linear", fill_value="extrapolate")(log.(Ep[Ep .>= E[end]]))]
    cross_section = pyconvert.(Array, cross_section)
    cross_section = vcat(cross_section[1], cross_section[2])
    cross_section = exp.(cross_section)
    cross_section[.!isfinite.(cross_section)] .= 0

    # Arbitrary correction of diss-cross-section
    correction_factor = 0.3 .+ 0.7 ./ (1 .+ exp.((Ep .- 120) ./ 20))
    cross_section = cross_section .* correction_factor

    cross_section[Ep .< 20.6] .= 0

    return cross_section
end

function e_N2ionx2sgp(Ep)
    # e_N2ionX2Sg+ - electron ionisation cross section (m^2) to the
    # ground-state  of N2+
    cross_section = similar(Ep)

    for iE in length(Ep):-1:1
        if Ep[iE] > 15.58 && Ep[iE] < 42.71
            cross_section[iE] = (1 - 15.58 / Ep[iE]) * exp(-834.7627 + 879.6264 * log(Ep[iE]) - 363.6978 * log(Ep[iE])^2 + 66.81782 * log(Ep[iE])^3 - 4.600032 * log(Ep[iE])^4)
        elseif Ep[iE] >= 42.71 && Ep[iE] < 300
            cross_section[iE] = exp(-100.7150 + 51.41718 * log(Ep[iE]) - 15.61180 * log(Ep[iE])^2 + 2.119856 * log(Ep[iE])^3 - 0.1089437 * log(Ep[iE])^4)
        elseif Ep[iE] >= 300
            cross_section[iE] = 9.119416e-15 * log(0.0275 * Ep[iE]) / Ep[iE]
        else
            cross_section[iE] = 0
        end
    end

    cross_section /= 1e4

    return cross_section
end

function e_N2iona2pu(Ep)
    # e_N2ionA2Pu - electron ionisation cross section (m^2) to the
    # first electronically excited state of N2+
    cross_section = similar(Ep)

    for iE in length(Ep):-1:1
        if Ep[iE] > 16.73 && Ep[iE] < 42.85
            cross_section[iE] = (1 - 16.73 / Ep[iE]) * exp(-40.81558 - 6.435371 * log(Ep[iE]) + 6.032905 * log(Ep[iE])^2 - 1.545984 * log(Ep[iE])^3 + 0.1261087 * log(Ep[iE])^4)
        elseif Ep[iE] >= 42.85 && Ep[iE] < 300
            cross_section[iE] = exp(-63.63995 + 19.39587 * log(Ep[iE]) - 5.333680 * log(Ep[iE])^2 + 0.6665464 * log(Ep[iE])^3 - 0.03254920 * log(Ep[iE])^4)
        elseif Ep[iE] >= 300
            cross_section[iE] = 9.022592e-15 * log(0.0275 * Ep[iE]) / Ep[iE]
        else
            cross_section[iE] = 0
        end
    end

    cross_section /= 1e4

    return cross_section
end

function e_N2ionb2sup(Ep)
    # e_N2ionB2Su+ - electron ionisation cross section (m^2) to the
    # second electronically excited state of N2+
    cross_section = similar(Ep)

    for iE in length(Ep):-1:1
        if Ep[iE] > 18.75 && Ep[iE] < 43.4
            cross_section[iE] = (1 - 18.75 / Ep[iE]) * exp(238.9714 - 290.0367 * log(Ep[iE]) + 112.9749 * log(Ep[iE])^2 - 19.41400 * log(Ep[iE])^3 + 1.241946 * log(Ep[iE])^4)
        elseif Ep[iE] >= 43.4 && Ep[iE] < 300
            cross_section[iE] = exp(-65.15851 + 19.39587 * log(Ep[iE]) - 5.333680 * log(Ep[iE])^2 + 0.6665464 * log(Ep[iE])^3 - 0.03254920 * log(Ep[iE])^4)
        elseif Ep[iE] >= 300
            cross_section[iE] = 1.976205e-15 * log(0.0275 * Ep[iE]) / Ep[iE]
        else
            cross_section[iE] = 0
        end
    end

    cross_section *= 1e-4

    return cross_section
end

function e_N2dion(Ep)
    # e_N2dion - dissociative ionization cross section (m^2)

    cross_section = similar(Ep)

    for iE in length(Ep):-1:1
        if Ep[iE] > 24.00 && Ep[iE] < 42.82
            cross_section[iE] = (1 - 24.00 / Ep[iE]) * exp(2066.673 - 2216.257 * log(Ep[iE]) + 870.8756 * log(Ep[iE])^2 - 151.3981 * log(Ep[iE])^3 + 9.830359 * log(Ep[iE])^4)
        elseif Ep[iE] >= 42.82 && Ep[iE] < 600
            cross_section[iE] = exp(-159.5311 + 90.55610 * log(Ep[iE]) - 25.06691 * log(Ep[iE])^2 + 3.076254 * log(Ep[iE])^3 - 0.1420719 * log(Ep[iE])^4)
        elseif Ep[iE] >= 600
            cross_section[iE] = 3.526978e-15 * log(0.0275 * Ep[iE]) / Ep[iE]
        else
            cross_section[iE] = 0
        end
    end

    cross_section *= 1e-4

    return cross_section
end

function e_N2ddion(Ep)
    # e_N2ddion - double dissociative ionization cross section (m^2)
    # i.e. N2 + e** -> e** + 2e* + N^+ + N^+

    cross_section = similar(Ep)

    for iE in length(Ep):-1:1
        if Ep[iE] > 42.00 && Ep[iE] < 44.45
            cross_section[iE] = (1 - 42.00 / Ep[iE]) * exp(830.9748 - 744.6075 * log(Ep[iE]) + 235.7621 * log(Ep[iE])^2 - 32.65055 * log(Ep[iE])^3 + 1.664205 * log(Ep[iE])^4)
        elseif Ep[iE] >= 44.45 && Ep[iE] < 550
            cross_section[iE] = exp(-240.8384 + 137.2195 * log(Ep[iE]) - 34.66862 * log(Ep[iE])^2 + 3.873074 * log(Ep[iE])^3 - 0.1624777 * log(Ep[iE])^4)
        elseif Ep[iE] >= 550
            cross_section[iE] = 9.873940e-16 * log(0.0275 * Ep[iE]) / Ep[iE]
        else
            cross_section[iE] = 0
        end
    end

    cross_section *= 1e-4

    return cross_section
end










# function e_N2d3sup(E)
#     cross_section = similar(E)

#     for ie in eachindex(E)
#         if E[ie] > 12.85 && E[ie] < 50
#             cross_section[ie] = (1 - 12.85 / E[ie]) * exp(0.938347 - 31.08899 * log(E[ie]) + 8.257418 * log(E[ie])^2 - 0.793736 * log(E[ie])^3)
#         elseif E[ie] >= 50
#             cross_section[ie] = 6.3166e-14 / E[ie]^3
#         else
#             cross_section[ie] = 0
#         end
#     end

#     return cross_section * 1e-4
# end
