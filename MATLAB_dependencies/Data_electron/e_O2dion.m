function cross_section = e_O2dion(Ep)
% Xs = e_O2dion(E)
% 
% e_O2dion - dissociative ionization cross section (m^2)
% E electron energy (eV)


for iE = numel(Ep):-1:1
  if Ep(iE) > 18.9 & Ep(iE) < 200,
    cross_section(iE) = (1-18.9/Ep(iE))*exp(-112.0394+59.82636*log(Ep(iE))-17.95648*log(Ep(iE))^2+2.411165*log(Ep(iE))^3-0.1228601*log(Ep(iE))^4);
  elseif Ep(iE) >= 200,
    cross_section(iE) = 6.078054e-15*log(0.030992*Ep(iE))/Ep(iE);
  else
    cross_section(iE) = 0;
  end
end
cross_section = cross_section/1e4;