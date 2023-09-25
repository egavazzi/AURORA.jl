function Xs = e_Ofine_2_0(E)
% Xs = e_Ofine_2_0(E)
% 
% e_Ofine_2_0 - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

xs = [2.1119   2.7641   6.779   8.8726   9.7054   9.7054   8.2582   7.7555 7.4153    6.779   6.4817   6.1974]*1e-22;%-18;
E_o = [0.0281 0.03 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.6 0.8 1];

Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
i = find(~isfinite(Xs));
Xs(i) = 0;
inans = find(E<0.0281);
Xs(inans) = 0;
