function Xs = e_N2vib0_3(Ep)
% Xs = e_N2vib0_3(E)
% 
% e_N2vib0_3- electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

E = [1.0271     1.8119     1.9281     2.3916     2.4163     2.5739 ...
     2.7995     3.0116     3.1689     3.3487     3.5064     3.7333 ...
     3.8524     4.0669     4.067      4.2041     4.3523     4.5029 ...
     4.6469]; 
s = [7.0427e-25 5.5034e-24 2.5652e-22 1.673e-20  1.6872e-20 3.4988e-21 ...
     1.5817e-20 2.1678e-21 1.1187e-20 1.5309e-21 4.3424e-21 1.2751e-21 ...
     1.7867e-21 7.6565e-22 7.6526e-22 1.0711e-21 5.0849e-22 1.0811e-21 ...
     5.5485e-22];
Xs = interp1(E,s,Ep,'pchip').*( Ep < E(end) );
Xs01 =  e_N2vib0_1(Ep).*(Ep>E(end));
Xs = Xs+Xs01;
inans = find(Ep<0.8559);
Xs(inans) = 0;
Xs(Ep>10) = 0; % TODO: FIX THIS/BG20190312