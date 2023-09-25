function Xs = e_O2_4p5(E)
% Xs = e_O2_4p5(E)
% 
% e_O2_4p5 - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

xs = [4.9802   6.8369   56.532   95.863   95.863   69.829   7.5986   1.7319]*1e-23;%-19;
E_o =[4.9      5        6         7        8       10      20       30];

Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
inans = find(~isfinite(Xs));
Xs(inans) = 0;
inans = find(E<4.9);
Xs(inans) = 0;

