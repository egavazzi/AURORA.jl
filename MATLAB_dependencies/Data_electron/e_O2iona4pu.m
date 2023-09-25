function cross_section = e_O2iona4pu(Energy)
% Xs = e_O2iona4pu(E)
% 
% e_O2iona4Pu - electron ionisation cross section (m^2) to the
% first excited state of O2+
% E electron energy (eV)


for iE = numel(Energy):-1:1
  if Energy(iE) > 16.1 & Energy(iE) < 300,
    cross_section(iE) = (1-16.1/Energy(iE))*exp(-49.02795+6.907601*log(Energy(iE))-1.325630*log(Energy(iE))^2+0.08157475*log(Energy(iE))^3-0.0002650063*log(Energy(iE))^4);
  elseif Energy(iE) >= 300,
    cross_section(iE) = 4.911588e-15*log(0.024825*Energy(iE))/Energy(iE);
  else
    cross_section(iE) = 0;
  end
end

cross_section = cross_section/1e4;
