function cross_section = e_N2ionb2sup(Energy)
% Xs = e_N2ionb2sup(E)
% 
% e_N2ionB2Su+ - electron ionisation cross section (m^2) to the
% second electronically excited state of N2+
% E electron energy (eV)


for iE = numel(Energy):-1:1,
  if Energy(iE) > 18.75 & Energy(iE) < 43.4,
    cross_section(iE) = (1-18.75/Energy(iE))*exp(238.9714-290.0367*log(Energy(iE))+112.9749*log(Energy(iE))^2-19.41400*log(Energy(iE))^3+1.241946*log(Energy(iE))^4);
  elseif Energy(iE) >= 43.4 & Energy(iE) < 300,
    cross_section(iE) = exp(-65.15851+19.39587*log(Energy(iE))-5.333680*log(Energy(iE))^2+0.6665464*log(Energy(iE))^3-0.03254920*log(Energy(iE))^4);
  elseif Energy(iE) >= 300,
    cross_section(iE) = 1.976205e-15*log(0.0275*Energy(iE))/Energy(iE);
  else
    cross_section(iE) = 0;
  end
end

cross_section = cross_section*1e-4;
