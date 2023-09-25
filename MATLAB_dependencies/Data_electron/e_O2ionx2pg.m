function cross_section = e_O2ionx2pg(Energy)
% Xs = e_O2ionx2pg(E)
% 
% e_O2ionX2Pg - electron ionisation cross section (m^2) to the
% ground-state of O2+
% E electron energy (eV)


for iE = numel(Energy):-1:1
  if Energy(iE) > 12.1 & Energy(iE) < 300,
    cross_section(iE) = (1-12.1/Energy(iE))*exp(-92.26781+44.59405*log(Energy(iE))-13.45450*log(Energy(iE))^2+1.809351*log(Energy(iE))^3-0.09207181*log(Energy(iE))^4);
  elseif Energy(iE) >= 300,
    cross_section(iE) = 8.771133e-15*log(0.024825*Energy(iE))/Energy(iE);
  else
    cross_section(iE) = 0;
  end
end

cross_section = cross_section/1e4;
