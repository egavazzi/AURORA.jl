function cross_section = e_N2cp3pu(Ep)
% Xs = e_N2c3pu(E)
% 
% e_N2c'3Pu - electron excitation cross section (m^2)
% E electron energy (eV)


for iE = numel(Ep):-1:1,
  if Ep(iE) > 12.08 & Ep(iE) < 18.36,
    cross_section(iE) = (1-12.08/Ep(iE))*exp(-3611.089+3701.604*log(Ep(iE))-1276.938*log(Ep(iE))^2+146.5513*log(Ep(iE))^3);
  elseif Ep(iE) >= 18.36 & Ep(iE) < 90,
    cross_section(iE) = exp(33.71613-61.7778*log(Ep(iE))+16.6116*log(Ep(iE))^2-1.50206*log(Ep(iE))^3);    
  elseif Ep(iE) >= 90,
    cross_section(iE) = 2.616e-14/Ep(iE)^3;
  else
    cross_section(iE) = 0;
  end
end

cross_section = cross_section*1e-4;
