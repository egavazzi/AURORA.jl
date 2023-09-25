function Xs = e_O2_8p4(E)
% Xs = e_O2_8p4(E)
% 
% e_O2_8p4 - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later


xs = [0.68369   9.5863    97.91   120.94   120.94    97.91   44.339    13.16]*1e-22;%-18;
E_o =[8.7992    9         10       20       50       90     120       600];

Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
inans = find(~isfinite(Xs));
Xs(inans) = 0;
inans = find(E<8.67);
Xs(inans) = 0;

