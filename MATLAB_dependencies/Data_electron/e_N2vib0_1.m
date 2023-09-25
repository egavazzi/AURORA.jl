function Xs = e_N2vib0_1(Ep)
% Xs = e_N2vib0_1(E)
% 
% e_N2vib0_1 - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright � Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

E = [0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
s = [1   3   3.8 4.2 4.6 5.5 6.7 8.1]*1e-23;%-19;
E = [E 1.6  1.7  1.75 1.8 1.85 1.9 1.95 2   2.05 2.1 2.15 2.2 2.23 2.3 2.35 2.4 2.45 2.5 2.55 2.6 2.65 2.73 2.8  2.85 2.9   3     3.05 3.075 3.14 3.2  3.3  3.4];
s = [s [0.28 0.51 1.2  1.9 2.6  4.5 5.7  4.2 2.2  1.4 2.5  5.0 5.8  4.4 3.6  1.9 1.5  2.5 4.5  3.6 1.67 1.21 2.78 3.20 2.04  0.88  1.35 1.92  1.58 0.79 0.97 0.56]*1e-20];%-16
E = [E  5   7.5 10  15   18   20   23   25   30   50   75];
s = [s  [6.5 3.1  1.4 4.1  7.5 19.4 12.1  7.2  2.4  1.4  0.67]*1e-22];%-18

Xs = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
I = find(~isfinite(Xs));
Xs(I) = 0;

Xs = Xs;

inans = find(Ep<0.2888);
Xs(inans) = 0;
Xs(Ep>10) = 0; % TODO: FIX THIS/BG20190312
