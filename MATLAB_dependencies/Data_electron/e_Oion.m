function Xs = e_Oion(E)
% Xs = e_Oion(E)
% 
% e_Oion - O electron ionization cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

xs =[ 1      3.45 6.9 11.3 13.8 11.9  9.76  8.41  6.9    5.12   2.32]*1e-21;%-17;
E_o = [13.618 20   30   50  100  200  300   400   600   1000   3000];

Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
i = find(~isfinite(Xs));
Xs(i) = 0;
i = find(E<13.618);
Xs(i) = 0;
