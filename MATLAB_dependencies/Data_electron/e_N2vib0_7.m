function Xs = e_N2vib0_7(Ep)
% Xs = e_N2vib0_7(E)
% 
%  e_N2vib0_7 - electron excitation cross section (m^2)
% E electron energy (eV)


%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

E = [2.0073       2.2342       2.3663       2.5946       2.8222       3.0483       3.2516       3.3851       3.5291       4.1979];
s = [6.317e-23   1.3348e-22   6.3896e-22   2.4423e-21   4.2059e-22   1.6918e-21   3.8452e-22   7.5752e-22   3.9773e-22   2.5062e-22];
Xs = interp1(E,s,Ep,'pchip').*( Ep < E(end) );
Xs01 =  e_N2vib0_1(Ep).*(Ep>E(end))/5;
Xs = Xs+Xs01;
inans = find(Ep<1.9);
Xs(inans) = 0;
Xs(Ep>10) = 0; % TODO: FIX THIS/BG20190312