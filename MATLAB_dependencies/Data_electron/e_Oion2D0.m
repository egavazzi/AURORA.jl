function Xs = e_Oion2D0(E)
% e_Oion - O electron ionization cross section to O^+(2P_0) (m^2)
% 
% Calling:
%  Xs = e_Oion2D0(E)
% Input:
%  E  - electron energy (eV), double array [1 x nE]
% Output:
%  Xs - collision cross section, double array [1 x nE]
% 
% Data source: Tima Sergienko, private communication

%  Copyright � Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later
Xs = 0*E;
for iE = numel(E):-1:1,
  if E(iE) > 16.9 & E(iE) <= 250,
    Xs(iE) = (1-16.9/E(iE))*exp(-38.80617-1.369768*log(E(iE))+1.103263*log(E(iE))^2-0.2226944*log(E(iE))^3+0.01331574*log(E(iE))^4);
  elseif E(iE) > 250,
    Xs(iE) = 4.109687e-15*log(0.032*E(iE))/E(iE);
  end
end

Xs = Xs*1e-4;
