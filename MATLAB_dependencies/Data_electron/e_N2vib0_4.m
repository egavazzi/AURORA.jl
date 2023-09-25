function Xs = e_N2vib0_4(Ep)
% Xs = e_N2vib0_4(E)
% 
% e_N2vib0_4 - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

Xs = (.3* e_N2vib0_1(Ep)+.35 * e_N2vib0_2(Ep))*2/5;
E = [1.6667     1.8867     2.2155     2.535      2.7699 ...
     3.049      3.1412     3.3004     3.4765     3.2972 ...
     3.4881     3.6861     3.8694     4.0781];
s = [1.9526e-23 5.5125e-23 1.1669e-20 7.5773e-22 1.0646e-20 ...
     4.0738e-21 4.0361e-21 2.293e-21  2.2845e-21 2.314e-21 ...
     2.2741e-21 1.2528e-21 9.8467e-22 5.0127e-22];
E = [1.1342 E];
s = [1e-27 s];
Xs = interp1(E,s,Ep,'pchip').*( Ep < E(end) );
Xs01 =  e_N2vib0_1(Ep).*(Ep>E(end));
Xs = Xs+0.35378*Xs01;

inans = find(Ep<1.1342);
Xs(inans) = 0;
Xs(Ep>10) = 0; % TODO: FIX THIS/BG20190312