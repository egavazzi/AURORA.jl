function Xs = e_O3sp3D0(E)
% Xs = e_O3sp3D0(E)
% 
% e_O3sp3D0 - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

xs = [  3  5.5 5   5.8 4.5 4   2.5]*1e-22;%-18;
E_o = [12.6 20  30  50 100 150 200];

Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
i = find(~isfinite(Xs));
Xs(i) = 0;
i = find(E<12.54);
Xs(i) = 0;
