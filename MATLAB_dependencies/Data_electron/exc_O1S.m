function X_o1s = exc_O1S(E)
% EXC_O1S - O1S excitation cross section
%   
% Calling:
%  Xs = exc_O1S(E)
% Input:
%  E - electron energy (eV), double array [1 x nE]
% Output:
%  Xs - excitation cross section (m^2), double array [1 x nE]

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

xs_o1s = [.1 .6 2.72 3.35 3.19 1.31]*1e-18*1e-4;
E_o1s = [4.19 5 7 10 20 30];
X_o1s = exp(interp1(log(E_o1s),log(xs_o1s),log(E),'linear','extrap'));
i = find(~isfinite(X_o1s));
X_o1s(E<4.19) = 0;
X_o1s(i) = 0;
