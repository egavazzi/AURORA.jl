function Xs = e_Oion4S0(E)
% e_Oion - O electron ionization cross section to O^+(4S_0) (m^2)
% 
% Calling:
%  Xs = e_Oion4S0(E)
% Input:
%  E  - electron energy (eV), double array [1 x nE]
% Output:
%  Xs - collision cross section, double array [1 x nE]
% 
% Data source: Tima Sergienko, private communication
% Xs = e_Oion(E)
% 
% e_Oion - O electron ionization cross section (m^2)
% E electron energy (eV)

%  Copyright � Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later
Xs = 0*E;

for iE = numel(E):-1:1,,
  if E(iE) > 13.6 & E(iE) <=  250,
    Xs(iE) = 0.35*(1-13.6/E(iE))*exp(-38.13225-1.957729*log(E(iE))+1.526543*log(E(iE))^2-0.3056663*log(E(iE))^3+0.01849928*log(E(iE))^4);
  elseif E(iE) > 250,
    Xs(iE) = 4.760656e-15*log(0.032*E(iE))/E(iE);
  end
end

Xs = Xs*1e-4;
