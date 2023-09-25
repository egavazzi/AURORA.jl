function cross_section = e_N2M1M2(Ep)
% Xs = e_N2M1M2(E)
% 
% e_N2M1(M2) - electron excitation cross section (m^2)
% E electron energy (eV)

for iE = numel(Ep):-1:1,
  if Ep(iE) > 13.15 & Ep(iE) < 24.,
    cross_section(iE) = (1-13.15/Ep(iE))*exp(115.6489-142.1746*log(Ep(iE))+44.0739*log(Ep(iE))^2-4.637555*log(Ep(iE))^3);
  elseif Ep(iE) >= 24. & Ep(iE) < 80,
    cross_section(iE) = exp(5.57238-37.398*log(Ep(iE))+10.292*log(Ep(iE))^2-0.97713*log(Ep(iE))^3);    
  elseif Ep(iE) >= 80,
    cross_section(iE) = 1.199e-13/Ep(iE)^3;
  else
    cross_section(iE) = 0;
  end
end

cross_section = cross_section*1e-4;
