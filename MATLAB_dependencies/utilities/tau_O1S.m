function tau = tau_O1S(nO,nO2,nN2,Tn,ne,Te)
% tau = tau_O1S(nO,nO2,nN2,Tn,ne,Te)
% 
% tau_O1D - effective O1D lifetime
% nO, nO2, nN2 - oxygen, molecular oxygen and molecular nitrogen
% concentration, Tn neutral temperature, ne electron concentration
% and Te electron temperature.
% 
% Einstein coefficients from Itikawa (1989), O2-quenching rate from
% Capetanakis (1993)


A_O1S = [2.732e-4,7.601e-2,1.215]; % Itikawa 1989
R_emi = 1/sum(A_O1S);
% Quenchin by N2
q_N2 = 0;

% O1S Quenching by O2, Capetanakis 1993
% q_O2 = 2.32e-2*exp((-812-1.82e-3*Tn.^2)./Tn);
R_gas = 8.3144621;
q_O2 = 6.2e-17*exp(-(5.65e3 - 4.18e-2*Tn.^2)./(R_gas*Tn));

Qn = (q_O2.*nO2)*1e-6;

tau = 1./( Qn + 1/R_emi );
