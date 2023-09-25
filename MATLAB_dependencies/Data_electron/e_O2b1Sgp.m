function Xs = e_O2b1Sgp(E)
% Xs = e_O2b1sgp(E)
% 
% e_O2b1sgp - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

xs = [1.9249 2.6425 9.3859 15.916 19.66 19.66 15.916  9.0932  5.535 1.0213  0.28755]*1e-23;%-19;
E_o =[1.9    2      3       4      5     7    10     20      30    90     120];

Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
inans = find(~isfinite(Xs));
Xs(inans) = 0;
inans = find(E<1.9);
Xs(inans) = 0;

