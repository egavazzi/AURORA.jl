function Xs = e_O3s3S0(E)
% Xs = e_O3s3S0(E)
% 
% e_O3s3S0 - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

xs =  [4    9.83 10.67  11.17 8.67 7.9 6.33 5.67 4.27]*1e-22;%-18;
E_o = [9.6 13.33 16.667 20   30   50 100  150  200];

Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
i = find(~isfinite(Xs));
Xs(i) = 0;
i = find(E<9.521);
Xs(i) = 0;
