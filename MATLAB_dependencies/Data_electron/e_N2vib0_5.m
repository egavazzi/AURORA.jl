function Xs = e_N2vib0_5(Ep)
% Xs = e_N2vib0_5(E)
% 
% e_N2vib0_5 - electron excitation cross section (m^2)
% E electron energy (eV)


%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later


E = [1.7459     2.0816     2.3099     2.4942     2.6043 ...
     2.8488     3.1042     3.2912     3.4976     3.7519 ...
     3.8671    ];
s = [6.8315e-24 2.1617e-23 2.9028e-21 3.9878e-21 6.2584e-21 ...
     2.6986e-22 3.9496e-21 2.6388e-22 1.613e-21  1.7546e-22 ...
     5.1621e-22];
s = [s(1)/3 s];
E = [1.4088 E];
Xs = interp1(E,s,Ep,'pchip').*( Ep < E(end) );
Xs01 =  e_N2vib0_1(Ep).*(Ep>E(end));
Xs = Xs+Xs01/10.663*2.5;

inans = find(Ep<1.4088);
Xs(inans) = 0;
Xs(Ep>10) = 0; % TODO: FIX THIS/BG20190312