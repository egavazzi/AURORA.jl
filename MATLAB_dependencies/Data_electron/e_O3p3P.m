function Xs = e_O3p3P(E)
% Xs = e_O3p3P(E)
% 
% e_O3p3P - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

xs = [5.1  7.8   4   2.9   1.1]*1e-22;%-18;
E_o = [13.5 20   30  50   100];

% $$$ Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
%Xs = exp(interp1(log(E_o),log(xs),log(E(E<=E_o(end))),'pchip'));
Xs = exp(interp1(log(E_o),log(xs),log(E(E<=60)),'pchip'));
%if any(E>E_o(end))
if any(E>60)
  
  %Xs = [Xs,xs(end)*E_o(end)/log(E_o(end))*log(E(E>E_o(end)))./E(E>E_o(end))];
  Xs = [Xs,e_O3p3P(59.9)*60/log(60)*log(E(E>60))./E(E>60)];
  
end
i = find(~isfinite(Xs));
Xs(i) = 0;
i = find(E<10.99);
Xs(i) = 0;
