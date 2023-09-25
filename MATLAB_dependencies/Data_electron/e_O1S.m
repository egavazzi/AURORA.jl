function Xs = e_O1S(E)
% Xs = e_O1S(E)
% 
% e_O1S - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

for iE = numel(E):-1:1,
  if E(iE) > 4.17 & E(iE) < 30,
    Xs(1,iE) = (1.-4.17/E(iE)).*exp(-77.00326+61.07961*log(E(iE))-37.16693*log(E(iE))^2+10.03347*log(E(iE))^3-1.021318*log(E(iE))^4);
  elseif E(iE) >= 30,
    Xs(1,iE) = 3.24e-14./E(iE).^3;
  end
end

Xs = Xs/1e4;
