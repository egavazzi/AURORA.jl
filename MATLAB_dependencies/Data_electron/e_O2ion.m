function Xs = e_O2ion(E)
% Xs = e_O2ion(E)
% 
% e_O2ion - O2 electron ionisation cross section (m^2)
% E electron energy (eV)


%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later


xs = [.41148   28.374   10.104   1.0422]*1e-21;%-17;
E_o =[12.708   97.938 1000   10000];

Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
inans = find(~isfinite(Xs));
Xs(inans) = 0;
inans = find(E<12);
Xs(inans) = 0;
