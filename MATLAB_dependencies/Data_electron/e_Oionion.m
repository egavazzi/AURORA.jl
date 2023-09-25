function Xs = e_Oionion(E)
% e_Oionion - O electron double ionization cross section (m^2)
% 
% Calling:
%  Xs = e_Oion(E)
% Input:
%  E  - electron energy (eV), double [1 x nE]
% Output:
%  Xs - cross-section for double-ionization
%
% Source of data: Tima Sergienko, private communication

%  Copyright � Bjorn Gustavsson 20181122, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later
Xs = 0*E;
for iE = numel(E):-1:1,
  if E(iE) > 28.5 & E(iE) <= 250,
    Xs(1,iE) = (1-28.5/E(iE))*exp(-8.790669-22.50029*log(E(iE))+6.626340*log(E(iE))^2-0.8663901*log(E(iE))^3+0.04147531*log(E(iE))^4);
  elseif E(iE) > 250,
    Xs(1,iE) = 2.471376e-15*log(0.032*E(iE))/E(iE);
  end
end

Xs = Xs*1e-4;