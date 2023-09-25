% DATA_ELECTRON V20180626
% Electron collision cross-sections, for major thermospheric constituents
% 
% 1, Molecular Nitrogen Collision cross sections
%
% 1.0 Elastic collision cross section
%
%   e_N2elastic              - Elastic collision cross-section
%
% 1.1 Rotation and Vibration excitation
%
%   e_N2rot0_2               - N_2 2nd rotational state
%   e_N2rot0_4               - N_2 4th rotational state
%   e_N2rot0_6               - N_2 6th rotational state
%   e_N2rot0_8               - N_2 8th rotational state
%   e_N2vib0_1               - N_2 1st vibrational state
%   e_N2vib0_2               - N_2 2nd vibrational state
%   e_N2vib0_3               - N_2 3rd vibrational state
%   e_N2vib0_4               - N_2 4th vibrational state
%   e_N2vib0_5               - N_2 5th vibrational state
%   e_N2vib0_6               - N_2 6th vibrational state
%   e_N2vib0_7               - N_2 7th vibrational state
%
% 1.2 Electronically excited states
%
%   e_N2a3sup                - N_2(A^3\sigma_u^+) cross section
%   e_N2b3pg                 - N_2(B^3\Pi_g) cross section
%   e_N2w3du                 - N_2(W^3\Delta_u) cross section
%   e_N2bp3sum               - N_2(B'^3\Sigma_u^-) cross section
%   e_N2ap1sum               - N_2(a'^1\Sigma_u^-) cross section
%   e_N2a1pg                 - N_2(a^1\Pi^+) cross section
%   e_N2w1du                 - N_2(w^1\Delta_u) cross section
%   e_N2c3pu                 - N_2(C^3\Pi_u) cross section
%   e_N2e3sgp                - N_2(E^3\Sigma_g^+) cross section
%   e_N2cp3pu                - N_2(c'^1\Pi_u) cross section
%   e_N2ab1sgp               - N_2(a''^1\Sigma_g^+) cross section
%   e_N2f3pu                 - N_2(F^3\Pi_u) cross section
%   e_N2g3pu                 - N_2(G^3\Delta_g) cross section
%   e_N2d3sup                - N_2(D^3\Sigma_u^+) cross section
%   e_N2bp1sup               - N_2(b'^1\Sigma_u^+) cross section
%   e_N2cp1sup               - N_2(c_4'^1\Sigma_u^+) cross section
%   e_N2o1pu                 - N_2(o_3^1\Pi_u)
%   e_N2M1M2                 - Higher states
%
% 1.3 Ionization and dissociation
%
%   e_N2ion                  - N_2 + e -> 2e + N_2^+
%   e_N2iona2pu              - N_2^+(A^2\Pi_u) 
%   e_N2ionb2sup             - N_2^+(B^2\Sigma_u^+)
%   e_N2ionx2sgp             - N_2^+(X^2\Sigma_u^+)
%   e_N2ddion                - N_2 + e -> 3e + N^+ + N^+
%   e_N2dion                 - N_2 + e -> 2e + N + N^+
%   e_N2dissociation         - N_2 + e -> N + N
% 
% 2, Molecular Oxygen Collision cross sections
%
% 2.0 Elastic collision cross section
%
%   e_O2elastic              - 
%
% 2.1 Rotation and Vibration excitation
%
%   e_O2vib                  - Sum of vibrational excitation cross-sections
%
% 2.2 Electronically excited states
%
%   e_O2a1Dg                 - O_2(a^1\Delta_g)
%   e_O2b1Sgp                - O_2(b^1\Sigma_g^+)
%   e_O2_4p5                 - O_2(~4.5  eV threshold)
%   e_O2_6                   - O_2(~6    eV threshold)
%   e_O2_8p4                 - O_2(~8.4  eV threshold)
%   e_O2_9p97                - O_2(~9.97 eV threshold)
%
% 2.3 Ionization and dissociation
%
%   e_O2ion                  - O_2 + e -> 2e + O_2^+
%   e_O2ionx2pg              - O_2^+(X^2\Pi_g)
%   e_O2iona4pu              - O_2^+(a^4\Pi_u)
%   e_O2ion16p9              - O_2^+(A^2\Pi_u, B^2\Sigma_g^-, ^2\Pi_u, c^4\Sigma_g^-)
%   e_O2ionb4sgm             - O_2^+(b^4\Sigma_g^-)
%   e_O2dion                 - O_2 + e -> 2e + O + O^+
%   e_O2ddion                - O_2 + e -> 2e + O^+ + O^+
%   e_O2_OO3S                - dissociation: O_2 + e -> e + O + O(3S)
% 
% 3, Atomic Oxygen Collision cross sections
%
% 3.0 Elastic collision cross section
%
%   e_Oelastic               - Elastic collision cross-section
%
% 3.1 Fine-structure collision cross sections
%
%   e_Ofine_1_0              - J=1 -> J=0
%   e_Ofine_2_0              - J=2 -> J=0
%   e_Ofine_2_1              - J=2 -> J=1
%
% 3.2 Electronically excited states
%
%   e_O1D                    - O(^1D)
%   e_O1S                    - O(^1S)
%   e_O3p3P                  - O(3p^3P)
%   e_O3p5P                  - O(3p^5P)
%   e_O3s3S0                 - O(3s^3S^0)
%   e_O3sp3D0                - O(3s'^3D^0)
%
% 3.3 Ionization
%
%   e_Oion                   - O + e -> 2e + O^+
% 
% 4 Secondary electron spectra
%
%   N2_e_2nd_dist            - secondary electron energy spectra ionisation of N2.
%   O2_e_2nd_dist            - secondary electron energy spectra for ionisation of O2. 
%   O_e_2nd_dist             - secondary electron energy spectra for ionisation of O.
%
% 5 Phase functions and Back-scattering ratios
%
%   phase_fcn_N2             - for N_2
%   phase_fcn_O              - for O
%   phase_fcn_O2             - for O_2
%   back_scattering_42stream - Back scatter ratios for electron 2-stream transport
%
% 6 Emission cross-sections
%   exc_4278                 - Xs = exc_4278(E)
%   exc_7774_O               - EI_7774_O - Xs = exc_7774_O(E)
%   exc_7774_O2              - s = exc_7774_O2(E)
%   exc_8446_O               - Xs = exc_8446(E)
%   exc_8446_O2              - Xs = exc_8446_O2(E)
%   exc_O1D                  - O1D excitation cross section
%   exc_O1S                  - O1S excitation cross section

%   Copyright © 2019 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
