function Xs = e_O3s5S0(E)
% Xs = e_O3s5S0(E)
% 
% e_O3s5S0 - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright � Tima Sergienko, Tima.Sergienko@irf.se
%  This is free software, licensed under GNU GPL version 2 or later

Xs = 0*E;

for iE = numel(E):-1:1,,
  if E(iE) > 9.14 & E(iE) < 23,
    Xs(1,iE) = (1.-9.14/E(iE)).*exp(-280.6036+256.7227*log(E(iE))-89.59541*log(E(iE))^2+10.22137*log(E(iE))^3);
  elseif E(iE) >= 23,
    Xs(1,iE) = 7.6e-15./E(iE).^3;
  end
end

Xs = Xs/1e4;
