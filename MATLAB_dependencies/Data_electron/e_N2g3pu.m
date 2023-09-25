function cross_section = e_N2g3pu(Ep)
% Xs = e_N2g3pu(E)
% 
% e_N2G3Pu - electron excitation cross section (m^2)
% E electron energy (eV)

for iE = numel(Ep):-1:1,
  if Ep(iE) > 12.8 & Ep(iE) < 18.84,
    cross_section(iE) = (1-12.8/Ep(iE))*exp(163.9467-189.6174*log(Ep(iE))+60.6054*log(Ep(iE))^2-6.615773*log(Ep(iE))^3);
  elseif Ep(iE) >= 18.84 & Ep(iE) < 80,
    cross_section(iE) = exp(6.881715-37.398*log(Ep(iE))+10.292*log(Ep(iE))^2-0.97713*log(Ep(iE))^3);    
  elseif Ep(iE) >= 80,
    cross_section(iE) = 4.441e-13/Ep(iE)^3;
  else
    cross_section(iE) = 0;
  end
end    

cross_section = cross_section*1e-4;
