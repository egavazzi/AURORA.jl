function Xs = e_N2vib0_6(Ep)
% Xs = e_N2vib0_6(E)
% 
% e_N2vib0_6 - electron excitation cross section (m^2)
% E electron energy (eV)


%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

E = [1.5176     2.143      2.2964     2.5275    2.8863     3.0358     3.1708     3.4811];
s = [1.0846e-25 9.2113e-24 4.9432e-23 4.459e-21 3.4105e-22 1.1193e-21 2.9207e-22 2.1537e-22];

Xs = interp1(E,s,Ep,'pchip').*( Ep < E(end) );
Xs01 =  e_N2vib0_1(Ep).*(Ep>E(end));
Xs = Xs+Xs01/20;
inans = find(Ep<1.68);
Xs(inans) = 0;
Xs(Ep>10) = 0; % TODO: FIX THIS/BG20190312