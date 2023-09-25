function Xs = e_Oion2P0(E)
% e_Oion - O electron ionization cross section to O^+(2P_0) (m^2)
% 
% Calling:
%  Xs = e_Oion2P0(E)
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
  if E(iE) > 18.6 & E(iE) <=  250,
    Xs(iE) = (1-18.6/E(iE))*exp(-34.62552-4.788984*log(E(iE))+2.031120*log(E(iE))^2-0.3347122*log(E(iE))^3+0.01838017*log(E(iE))^4);
  elseif E(iE) > 250,
    Xs(iE) = 2.330830e-15*log(0.032*E(iE))/E(iE);
  end
end

Xs = Xs*1e-4;
