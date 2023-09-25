function Xs = e_O2vib(E)
% Xs = e_O2vib(E)
% 
% e_O2vib - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

xs = [1.9973   17.464   26.395   23.806   1.255  2.5856   11.555   1.6247]*1e-22;%-18;
E_o =[0.24755      0.44367       0.6389      0.73923       1.4251       6.1282        9.086       14.7];

Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
inans = find(~isfinite(Xs));
Xs(inans) = 0;
inans = find(E<0.24);
Xs(inans) = 0;
