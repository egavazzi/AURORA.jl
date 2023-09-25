function Xs = e_N2vib0_2(Ep)
% Xs = e_N2vib0_2(E)
% 
% e_N2vib0_2 - electron excitation cross section (m^2)
% E electron energy (eV)

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

E = [0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
s = [1   3   3.8 4.2 4.6 5.5 6.7 8.1]*1e-23;%-19;
E = [E 1.67       1.75       1.8       1.85       1.9       1.9375        1.978       2.0  2.05        2.1       2.125       2.2       2.24       2.28       2.35       2.4 2.45   2.5     2.55       2.6       2.65       2.7       2.75       2.8       2.85 2.9         2.95        3         3.05       3.1       3.15       3.2       3.25 3.3       3.35         3.4];
s = [s [ 0.08      0.19367      0.28       1.37       1.94       3.82       4.04       3.65 3.16       2.24       1.79      0.57       1.14       2.47       4.16      3.88 2.16       1.25      0.51      0.97       1.94       2.32       1.67      0.97 0.4      0.51      0.91       1.2       1.08      0.39      0.23      0.54 0.74       0.46      0.17      0.20]*1e-20];%-16
E = [E   5   7.5 10  15   18   20   23   25   30   50   75];
s = [s  [6.5 3.1  1.4 4.1  7.5 19.4 12.1  7.2  2.4  1.4  0.67]*1e-22];%-18

Xs = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
I = find(~isfinite(Xs));
Xs(I) = 0;

Xs = Xs;

inans = find(Ep<0.2888);
Xs(inans) = 0;
Xs(Ep>10) = 0; % TODO: FIX THIS/BG20190312
