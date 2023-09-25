function Xs = e_O2_9p97(E)
% Xs = e_O2_9p97(E)
% 
% e_O2_9p97 - electron excitation cross section (m^2)
% E electron energy (eV)


%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later


xs = [0.12616   1.4321   3.7052   5.9598    6.9829   6.283   3.5146    2.185    1.316]*1e-22;%-18;
E_o =[10.66    20       40       60        90      100     200       400      600];

Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
inans = find(~isfinite(Xs));
Xs(inans) = 0;
inans = find(E<9.97);
Xs(inans) = 0;

