function Xs = e_O2_OO3S(E)
% Xs = e_O2_OO3S(E)
% 
% e_O2_OO3S - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

xs = [ 0.15583   1.2885    2.185   2.9369   3.3692   3.3338   2.6989   1.7689]*1e-22;%-18;
E_o =[15.647    18.361    30      50       80      100      200      400];

Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
inans = find(~isfinite(Xs));
Xs(inans) = 0;
inans = find(E<15.6);
Xs(inans) = 0;

