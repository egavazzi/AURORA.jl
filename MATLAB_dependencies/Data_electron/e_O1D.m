function Xs = e_O1D(E)
% Xs = e_O1D(E)
% 
% e_O1D - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

for iE = numel(E):-1:1,
  if E(iE) > 1.967 & E(iE) < 6.867,
    Xs(1,iE) = (1.-1.9./E(iE)).*exp(-38.0685-0.2992.*E(iE)+0.20375*E(iE).^2-0.0211739.*E(iE).^3);
  elseif E(iE) > 6.867 & E(iE) < 30,
    Xs(1,iE) = exp(-34.081-0.912397.*E(iE)+7.185417e-2.*E(iE).^2-2.48398e-3.*E(iE).^3+3.00574e-5.*E(iE).^4);
  elseif E(iE) >=  30,
    Xs(1,iE) = 1.881e-13./E(iE).^3;
  end
end
Xs = Xs/1e4;