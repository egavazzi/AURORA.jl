function Xs = e_O2a1Dg(E)
% Xs = e_O2a1dg(E)
% 
% e_O2a1dg - electron excitation cross section (m^2)
% E electron energy (eV)


%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later


xs = [6.1516   15.916   37.052   56.532   73.616   86.254   89.031   75.986 33.338    19.66    5.535   2.9369]*1e-23;%-19;
E_o =[1.4678    2        3        4        5        6        7        9     20        30      80     100];

Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
inans = find(~isfinite(Xs));
Xs(inans) = 0;
inans = find(E<0.977);
Xs(inans) = 0;
