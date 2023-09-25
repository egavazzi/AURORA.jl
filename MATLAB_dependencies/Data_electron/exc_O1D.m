function X_o1d = exc_O1D(E)
% EXC_O1D - O1D excitation cross section
%   
% Calling:
%  Xs = exc_O1D(E)
% Input:
%  E - electron energy (eV), double array [1 x nE]
% Output:
%  Xs - excitation cross section (m^2), double array [1 x nE]

%  Copyright © Bjorn Gustavsson 20110527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

xs_o1d = [3.02 4.42 5.43 3.89 3.44 2.57 2.02 2.72 2.07 1.99 0.8 0.7]*1e-17;
E_o1d = [4 5 6 7 7.1 9 9.6 10 15.6 20 28.1 30];
xs_o1di = [0.35 1.6 3.02 4.42 5.43 mean([3.89 3.44]) mean([2.57 2.02 2.72]) 1.8 1.5 0.8 0.7]*1e-17;
E_o1di =  [2     3   4    5    6  7.05                9.6                   15.6 20 28.1 30];
X_o1d = exp(interp1(log(E_o1di),log(xs_o1di),log(E),'linear','extrap'))*1e-4;
inans = find(~isfinite(X_o1d));
X_o1d(inans) = 0;
