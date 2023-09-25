function Xs = e_O2elastic(E)
% Xs =  e_O2elastic(E)
% 
% e_O2elastic - elastic electron collision cross section (m^2)
% E electron energy (eV)


%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later


xs =[ 4.0725   4.2883    5.272   5.8454   7.1862   8.8346   9.7956   10.861   8.8346   6.4813   4.5155   2.2608    0.91126  0.10422]*1e-20;
E_o =[0.05     0.1       0.6     1        2        7       10        13.58   30       60      100      300      1000        1e4];
Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
i = find(~isfinite(Xs));
Xs(i) = 0;
i = find(E<1e-2);
Xs(i) = 0;
