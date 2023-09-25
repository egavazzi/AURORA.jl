function cross_section = e_N2f3pu(Ep)
% Xs = e_N2F3pu(E)
% 
% e_N2F3Pu - electron excitation cross section (m^2)
% E electron energy (eV)

for iE = numel(Ep):-1:1,
  if Ep(iE) > 12.75 & Ep(iE) < 23.28,
    cross_section(iE) = (1-12.75/Ep(iE))*exp(27.4082-52.0658*log(Ep(iE))+13.78086*log(Ep(iE))^2-1.26079*log(Ep(iE))^3);
  elseif Ep(iE) >=  23.28 & Ep(iE) < 90,
    cross_section(iE) = exp(1.13172-32.235*log(Ep(iE))+8.5621*log(Ep(iE))^2-0.78723*log(Ep(iE))^3);    
  elseif Ep(iE) >= 90,
    cross_section(iE) = 3.171e-13/Ep(iE)^3;
  else
    cross_section(iE) = 0;
  end
end

cross_section = cross_section*1e-4;
