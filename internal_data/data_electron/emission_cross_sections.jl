using DataInterpolations: LinearInterpolation, PCHIPInterpolation
using Statistics: mean

"""
    excitation_4278(E)

Returns the electron-impact excitation cross-section for 4278Å from N2 ions.

Best fit to experiments (Tima Sergienko, personal communication).

# Calling
σ = excitation\\_4278(E)

# Input
- `E`: incoming electron energy (eV)

# Output
- `σ`: emission cross-section (m²)
"""
function excitation_4278(E)
    σ = similar(E)
    for i_E in eachindex(E)
        if 18.75 < E[i_E] < 43.4
            σ[i_E] = (1 - 18.75 / E[i_E]) *
                                 exp(238.9714 - 290.0367 * log(E[i_E]) +
                                     112.9749 * log(E[i_E])^2 - 19.41400 * log(E[i_E])^3 +
                                     1.241946 * log(E[i_E])^4)
        elseif 43.4 <= E[i_E] < 300
            σ[i_E] = exp(-65.15851 + 19.39587 * log(E[i_E]) -
                                     5.333680 * log(E[i_E])^2 + 0.6665464 * log(E[i_E])^3 -
                                     0.03254920 * log(E[i_E])^4)
        elseif E[i_E] >= 300
            σ[i_E] = 1.976205e-15 * log(0.0275 * E[i_E]) / E[i_E]
        else
            σ[i_E] = 0
        end
    end

    # Scale the result
    σ = σ / 1e4 / 3.3773

    return σ
end

"""
    excitation_6730_N2(E)

Returns the electron-impact excitation cross-section for 6730Å bands from N2 for transitions 4-1 and 5-2.

Digitised and extrapolated from Lanchester et al. (2009), p. 2545.
https://doi.org/10.5194/angeo-27-2881-2009

% parent       N2
% products     0    (emission)
% threshold    7.3532   (???)
% units        eV

# Calling
σ = excitation\\_6730\\_N2(E)

# Input
- `E`: incoming electron energy (eV)

# Output
- `σ`: emission cross-section (m²)
"""
function excitation_6730_N2(E)
    # Data points for interpolation (first column is energy and second column is corresponding cross-section)
    EnX = [
        6.37543473010668    7.96983399250628e-24;
        6.645257176389      1.06190815625254e-22;
        7.05092034556325    1.64964807409802e-22;
        7.4578382213999     2.16886182213227e-22;
        8.00449957408126    2.75790759057769e-22;
        8.90007900596942    3.79935736503941e-22;
        9.49542375014288    4.40025921260209e-22;
        9.95343543588585    4.9619476030029e-22;
        10.5256226847025    5.52114287040292e-22;
        10.8077889919725    5.70849650019635e-22;
        11.064999704106     5.90220775981765e-22;
        11.4951949205798    5.8629466459164e-22;
        12.0474108103454    5.70849650019635e-22;
        12.6997233300692    5.30441727984392e-22;
        13.2701234933231    4.9619476030029e-22;
        14.5729827646183    4.06158598837698e-22;
        16.5281054399492    3.2154841769777e-22;
        19.4730095831928    2.39723243377653e-22;
        25.0476823843719    1.46290926796731e-22;
        28.406890677532     1.12014788271289e-22;
        31.1047253633823    9.23033346932114e-23;
        33.6615955341341    7.65698808012419e-23;
        36.0050530813703    6.56736727188188e-23;
        39.3076282437633    5.2691326305265e-23;
        42.6638138766614    4.37098896574342e-23;
        48.1015886772302    3.30248291643171e-23;
        53.759362034737     2.57984849916881e-23;
        91.0396256375159    8.18546730706903e-24
    ]

    # Convert data to log-log space
    log_energy = log.(EnX[:, 1])
    log_cross_section = log.(EnX[:, 2])

    # Create the interpolation function
    interp = LinearInterpolation(log_cross_section, log_energy; extrapolation = ExtrapolationType.Extension)

    # Compute interpolated values
    σ = exp.(interp(log.(E)))

    # Set to 0 under threshold
    E_th = 7.3532
    σ[E .< E_th] .= 0.0

    return σ
end

"""
    excitation_8446_O(E)

Returns the electron-impact excitation cross-section for 8446Å from O.
Produced by the transition
    ground-state --(different pathways)--> OI 3p 3P --(prompt emission of 8446Å)--> ??

From Itikawa and Ichimura (1990).
https://doi.org/10.1063/1.555857

% parent       O
% products     0    (emission)
% threshold    10.99   (Itikawa) Confirmed/BG-20191016
% units        eV

# Calling
σ = excitation\\_8446\\_O(E)

# Input
- `E`: incoming electron energy (eV)

# Output
- `σ`: emission cross-section (m²)
"""
function excitation_8446_O(E)
    # Data points for interpolation (first column is energy and second column is corresponding cross-section)
    EnX = [
        1.1000000e+01 2.5781459e-01;
        1.2059649e+01 3.6743271e-01;
        1.3221375e+01 4.8399335e-01;
        1.4495013e+01 5.9391952e-01;
        1.5891342e+01 6.8434932e-01;
        1.7422182e+01 7.4632073e-01;
        1.9100490e+01 7.7643659e-01;
        2.0940473e+01 7.6493066e-01;
        2.2957704e+01 6.7097864e-01;
        2.5169259e+01 5.5059069e-01;
        2.7593856e+01 4.5203505e-01;
        3.0252019e+01 3.9695472e-01;
        3.3166248e+01 3.6965228e-01;
        3.6361209e+01 3.5009898e-01;
        3.9863946e+01 3.3411540e-01;
        4.3704108e+01 3.1832706e-01;
        4.7914199e+01 2.9997511e-01;
        5.2529855e+01 2.7743037e-01;
        5.7590145e+01 2.5248291e-01;
        6.3137901e+01 2.3229768e-01;
        6.9220083e+01 2.1658731e-01;
        7.5888171e+01 2.0184430e-01;
        8.3198607e+01 1.8801995e-01;
        9.1213271e+01 1.7506665e-01;
        1.0000000e+02 1.6293805e-01;
        2.1544347e+02 8.8234000e-02;
        4.6415888e+02 4.6805251e-02;
        1.0000000e+03 2.4440707e-02;
        2.1544347e+03 1.2604857e-02;
        4.6415888e+03 6.4357220e-03;
        1.0000000e+04 3.2587610e-03;
        2.1544347e+04 1.6386314e-03;
        4.6415888e+04 8.1909190e-04;
        1.0000000e+05 4.0734512e-04
    ]

    # Convert data to log-log space
    log_energy = log.(EnX[:, 1])
    log_cross_section = log.(EnX[:, 2])

    # Create the interpolation function
    interp = PCHIPInterpolation(log_cross_section, log_energy; extrapolation = ExtrapolationType.Extension)

    # Compute interpolated values
    σ = exp.(interp(log.(E)))

    # Set to 0 under threshold
    E_th = 10.99
    σ[E .< E_th] .= 0.0

    # Convert units
    σ .*= 1e-21

    return σ
end

"""
    excitation_8446_O2(E)

Returns the electron-impact excitation cross-section for 8446Å from O2.
Produced by the dissociative ionization/excitation of
    O2 --(different pathways)--> OI 3p 3P --(prompt emission of 8446Å)--> ??

From Schulman et al. (1985).
Digitised from fig.7 by Daniel Whiter in March 2009. Value in paper is 2 +/- 15% @ 100eV.
https://doi.org/10.1103/PhysRevA.32.2100

% parent       O2
% products     0    (emission)
% threshold    15.9   (???) - seems low, dissociation energy of O2
%                             is 5.16 eV and the energy-level of
%                             O(3p3P) is 10.99 eV, so threshold has
%                             to  be 16.15 eV /BG 20180529

# Calling
σ = excitation\\_8446\\_O2(E)

# Input
- `E`: incoming electron energy (eV)

# Output
- `σ`: emission cross-section (m²)
"""
function excitation_8446_O2(E)
    # Data points for interpolation (first column is energy and second column is corresponding cross-section)
    EnX = [
        15.9877      0.00888889;
        16.3801      0.124444;
        16.7724      0.24;
        17.1341      0.377778;
        18.1027      0.475556;
        19.0345      0.6;
        21.5111      0.804444;
        24.0429      0.968889;
        26.5992      1.11556;
        31.3502      1.27111;
        34.5563      1.34667;
        37.787       1.40444;
        41.5816      1.45333;
        42.6299      1.49333;
        44.7387      1.56444;
        48.4843      1.64889;
        51.1387      1.72444;
        57.6123      1.83111;
        64.6682      1.91556;
        70.6697      1.96444;
        77.2536      1.99111;
        84.9594      2.00444;
        93.2475      1.99556;
        103.767      1.96889;
        114.863      1.92444;
        130.973      1.84444;
        147.09       1.76;
        159.877      1.68889;
        178.207      1.6;
        195.96       1.52889;
        215.375      1.45333;
        231.448      1.4;
        253.051      1.33778;
        279.08       1.26667;
        301.769      1.21778;
        327.221      1.16444;
        348.794      1.12444;
        375.344      1.07556;
        398.572      1.03556;
        1e3          0.58735;
        1e4          0.14193;
        1e5          0.034298;
    ]
    # The last 3 rows were added by B. Gustavsson 20180625 to make extrapolation to higher energies well-behaved.
    # These values are based on polynomial extrapolation of log of cross-section vs log of energy

    # Convert data to log-log space
    log_energy = log.(EnX[:, 1])
    log_cross_section = log.(EnX[:, 2])

    # Create the interpolation function
    interp = PCHIPInterpolation(log_cross_section, log_energy; extrapolation = ExtrapolationType.Extension)

    # Compute interpolated values
    σ = exp.(interp(log.(E)))

    # Set to 0 under threshold
    E_th = 16.15
    σ[E .< E_th] .= 0.0

    # Convert units
    σ .*= 1e-22

    return σ
end

"""
    excitation_7774_O(E)

Returns the electron-impact excitation cross-section for 7774Å from O.
Produced by the transition
    groud-state --(different pathways)--> OI 3p 5P --(prompt emission of 7774Å)--> OI 3s 5S

Digitized and extrapolated from Julienne and Davis (1976), p. 1397.
https://doi.org/10.1029/JA081i007p01397

% 	typed in by Mina Ashrafi July 2008
% parent       O
% products     0    (emission)
% threshold    10.74   (Itikawa Corrected:BG-20191016)
% units        eV

# Calling
σ = excitation\\_7774\\_O(E)

# Input
- `E`: incoming electron energy (eV)

# Output
- `σ`: emission cross-section (m²)
"""
function excitation_7774_O(E)
    # Data points for interpolation (first column is energy and second column is corresponding cross-section)
    EnX = [
        10.74      0.1;
        11.0929    0.7692;
        11.6841    1.3736;
        12.2774    2.5162;
        12.4855    3.5089;
        12.6895    4.6491;
        12.8848    5.5134;
        13.0787    6.4279;
        13.4525    7.3044;
        13.8208    7.7530;
        14.1878    8.0900;
        14.7353    8.2986;
        15.4639    8.4398;
        16.3723    8.3663;
        16.9164    8.2237;
        17.4598    8.0150;
        18.0032    7.8115;
        18.5466    7.6132;
        19.0901    7.4199;
        19.6328    7.1701;
        19.9937    6.9290;
        20.5372    6.7531;
        20.8981    6.5260;
        21.2590    6.3066;
        21.8025    6.1465;
        21.9816    5.9401;
        22.5243    5.7401;
        23.0664    5.4998;
        23.9701    5.1359;
        24.5122    4.9209;
        25.0543    4.7149;
        25.4145    4.5176;
        25.9566    4.3285;
        26.4987    4.1472;
        27.2212    3.9062;
        27.9438    3.6792;
        28.4852    3.4952;
        29.0272    3.3489;
        29.5680    3.1544;
        30.2905    2.9711;
        30.6515    2.8712;
        31.3747    2.7275;
        32.0972    2.5690;
        32.6386    2.4405;
        33.3618    2.3184;
        33.9032    2.2024;
        34.4453    2.1102;
        35.1685    2.0046;
        35.8911    1.8881;
        36.6136    1.7784;
        37.3362    1.6750;
        38.0594    1.5912;
        38.9631    1.4859;
        39.8668    1.3876;
        40.9516    1.2848;
        41.8560    1.2100;
        42.7603    1.1397;
        43.4829    1.0734;
        44.3872    1.0110;
        45.2923    0.9604;
        46.1959    0.8968;
        47.1003    0.8447;
        48.0047    0.7956;
        48.9097    0.7557;
    	10000	   0.1e-4; # Dummy-value
        ]

    # Convert data to log-log space
    log_energy = log.(EnX[:, 1])
    log_cross_section = log.(EnX[:, 2])

    # Create the interpolation function
    interp = PCHIPInterpolation(log_cross_section, log_energy; extrapolation = ExtrapolationType.Extension)

    # Compute interpolated values
    σ = exp.(interp(log.(E)))

    # Set to 0 under threshold
    E_th = 10.74
    σ[E .< E_th] .= 0.0

    # Convert units
    σ .*= 1e-22

    return σ
end

"""
    excitation_7774_O2(E)

Returns the electron-impact excitation cross-section for 7774Å from O2.
Produced by the dissociative ionization/excitation of
    O2 --(different pathways)--> OI 3p 5P --(prompt emission of 7774Å)--> OI 3s 5S

From Erdman and Zipf (1987), p. 4540.
https://doi.org/10.1063/1.453696

% 	typed in by Mykola/Nickolay Ivshenko August 2006
% parent       O2
% products     0    (emission)
% threshold    15.9   (???) Seems OK, 10.74 + 5.15 (O2-bond-energy)
% units        eV

# Calling
σ = excitation\\_7774\\_O2(E)

# Input
- `E`: incoming electron energy (eV)

# Output
- `σ`: emission cross-section (m²)
"""
function excitation_7774_O2(E)
    # Data points for interpolation (first column is energy and second column is corresponding cross-section)
    EnX =  [
    16	    0.044;
	17	    0.35;
	18	    0.618;
	19	    0.888;
	20	    1.10;
	22	    1.52;
	24	    1.75;
	26	    1.91;
	28	    2.05;
	30	    2.16;
	35	    2.40;
	40	    2.62;
	50	    3.19;
	60	    3.83;
	70	    3.99;
	80	    4.07;
	90	    4.31;
	100	    4.23;
	125	    3.87;
	150	    3.42;
	175	    3.06;
	200	    2.76;
	225	    2.56;
	250	    2.34;
	275	    2.20;
	300	    2.04;
	325	    1.94;
	350	    1.82;
	375	    1.72;
	400	    1.65;
	425	    1.58;
	450	    1.50;
	475	    1.44;
	500	    1.39;
	600	    1.20;
	700	    1.06;
	800	    0.954;
	900	    0.868;
	1000	0.797;
	1200	0.687;
	1400	0.605;
	1600	0.542;
	1800	0.492;
	2000	0.451;
	2200	0.416;
	2400	0.387;
	2600	0.362;
	2800	0.340;
	3000	0.321;
	4000	0.251;
	5000	0.208;
	6000	0.178;
	7000	0.156;
	8000	0.139;
	9000	0.125;
	10000	0.114;
    ]

    # Convert data to log-log space
    log_energy = log.(EnX[:, 1])
    log_cross_section = log.(EnX[:, 2])

    # Create the interpolation function
    interp = PCHIPInterpolation(log_cross_section, log_energy; extrapolation = ExtrapolationType.Extension)

    # Compute interpolated values
    σ = exp.(interp(log.(E)))

    # Set to 0 under threshold
    E_th = 16.14
    σ[E .< E_th] .= 0.0

    # Convert units
    σ .*= 1e-22

    return σ
end

"""
    excitation_O1D(E)

Returns the electron-impact excitation cross-section for O1D (will emit 6300Å, not prompt).

From Itikawa and Ichimura (1990).
https://doi.org/10.1063/1.555857

# Calling
σ = excitation\\_O1D(E)

# Input
- `E`: incoming electron energy (eV)

# Output
- `σ`: emission cross-section (m²)
"""
function excitation_O1D(E)
    # Data points for interpolation
    xs_o1d = [0.35, 1.6, 3.02, 4.42, 5.43, mean([3.89, 3.44]), mean([2.57, 2.02, 2.72]), 1.8, 1.5, 0.8, 0.7] .* 1e-17 .* 1e-4
    E_o1d =  [2, 3, 4, 5, 6, 7.05, 9.6, 15.6, 20, 28.1, 30]

    # Convert data to log-log space
    log_energy = log.(E_o1d)
    log_cross_section = log.(xs_o1d)

    # Create the interpolation function
    interp = LinearInterpolation(log_cross_section, log_energy; extrapolation = ExtrapolationType.Extension)

    # Compute interpolated values
    σ = exp.(interp(log.(E)))

    # Set nonfinite values to 0
    I = findall(.!isfinite.(σ))
    σ[I] .= 0

    # Set to 0 under threshold
    E_th = 1.967
    σ[E .< E_th] .= 0.0

    return σ
end

"""
    excitation_O1S(E)

Returns the electron-impact excitation cross-section for O1S (will emit 5577Å, not prompt).

From Itikawa and Ichimura (1990).
https://doi.org/10.1063/1.555857

# Calling
σ = excitation\\_O1S(E)

# Input
- `E`: incoming electron energy (eV)

# Output
- `σ`: emission cross-section (m²)
"""
function excitation_O1S(E)
    # Data points for interpolation
    xs_o1s = [.1, .6, 2.72, 3.35, 3.19, 1.31] .* 1e-18 .* 1e-4
    E_o1s = [4.19, 5, 7, 10, 20, 30]

    # Convert data to log-log space
    log_energy = log.(E_o1s)
    log_cross_section = log.(xs_o1s)

    # Create the interpolation function
    interp = LinearInterpolation(log_cross_section, log_energy; extrapolation = ExtrapolationType.Extension)

    # Compute interpolated values
    σ = exp.(interp(log.(E)))

    # Set nonfinite values to 0
    I = findall(.!isfinite.(σ))
    σ[I] .= 0

    # Set to 0 under threshold
    E_th = 4.19
    σ[E .< E_th] .= 0.0

    return σ
end
