function Xs = e_O3p5P(E)
% Xs = e_O3p5P(E)
% 
% e_O3p5P - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

% $$$ xs =[ 0.235  2.35      3.5    2.06 0.69 0.31 0.125]*1e-22;%-18;
% $$$ E_o=[10.74  15    15.625 20   30   50   70];
% $$$ 
% $$$ Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
% $$$ i = find(~isfinite(Xs));
% $$$ Xs(i) = 0;
% $$$ i = find(E<10.74);
% $$$ Xs(i) = 0;
Xs = 0*E;

for iE = numel(E):-1:1,,
  if E(iE) > 10.73 & E(iE) < 31.614,
    Xs(1,iE) = (1.-10.73/E(iE)).*exp(678.5487-931.9045*log(E(iE))+452.4033*log(E(iE))^2-97.15600*log(E(iE))^3+7.763117*log(E(iE))^4);
  elseif E(iE) > 31.614 & E(iE) < 200,
    Xs(1,iE) = exp(-54.32016+19.35011*log(E(iE))-8.191310*log(E(iE))^2+1.323729*log(E(iE))^3-0.07969778*log(E(iE))^4);
  elseif E(iE) >= 200,
    Xs(1,iE) = 1.593536e-14/E(iE)^3;
  end
end
Xs = Xs/1e4;
