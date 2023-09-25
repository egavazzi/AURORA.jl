function Xs = e_O2_6(E)
% Xs = e_O2_6(E)
% 
% e_O2_6 - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

xs = [59.598   146.26   235.27   223.16   200.79   56.532   45.767   28.453    19.66   9.3859   6.8369]*1e-23;%-19;
E_o =[6 7 8 9 10 17.783 20 30 40 90 100];

Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
inans = find(~isfinite(Xs));
Xs(inans) = 0;
inans = find(E<6);
Xs(inans) = 0;

