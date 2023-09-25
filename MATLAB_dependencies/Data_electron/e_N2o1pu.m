function cross_section = e_N2o1pu(Ep)
% Xs = e_N2o1pu(E)
% 
% e_N2o1Pu - electron excitation cross section (m^2)
% E electron energy (eV)

for iE = numel(Ep):-1:1,
  if Ep(iE) > 13.1 & Ep(iE) < 150,
    cross_section(iE) = (1-13.1/Ep(iE))*exp(-44.40084+3.656131*log(Ep(iE))-0.9243795*log(Ep(iE))^2+0.0642402*log(Ep(iE))^3);
  elseif Ep(iE) >= 150,
    cross_section(iE) = 8.0809e-17*log(0.0573*Ep(iE))/Ep(iE);
  else
    cross_section(iE) = 0;
  end
end

cross_section = cross_section*1e-4;
