function tau = tau_O1D(nO,nO2,nN2,ne,Tn,Te,h)
% TAU_O1D - function for calulating O(1D) lifetime
% 
% Calling:
%  tau = tau_O1D(nO,nO2,nN2,ne,Tn,Te,h)
% 
% References for Quenching rates (Solonon et al 1988): 
%  N2 + O1D -> O3P + N2 [Streit et al 1976]
%  O2 + O1D -> O3P + N2 [Streit et al 1976]
%  O  + O1D -> O3P + N2 [Abreu 1986]
%  e- + O1D -> O3P + e- [Berrington and Burke, 1981]


% Quenching by N2
q7 = 2.0e-11*exp(107.8./Tn);

% Quenching by O2
q8 = 2.9e-11*exp(67.5./Tn);

% Quenching by thermal electrons
q9 = 1.6e-12*Te.^.91;

Qe = q9.*ne;

% Quenching by O
q10 = 8e-12;
% kolla med Tima!

% For calculation of exitation/emission/quenching coefficients
Qn =    q7.*nN2 + q8.*nO2 + q10*nO;

R = 1/107.5;

tau = 1./(Qn+Qe+R);
