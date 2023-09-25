function Xs = e_Ofine_1_0(E)
% Xs = e_Ofine_1_0(E)
% 
% e_Ofine_1_0 - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

xs = [ 0.61974   1.6626    2.176   2.4233   2.3803   2.2758    2.176   1.9893]*1e-21;%-17;
E_o = [0.02 0.03 0.04 0.06 0.08 0.1 0.3 100];

Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
i = find(~isfinite(Xs));
Xs(i) = 0;
inans = find(E<0.0281);
Xs(inans) = 0;
