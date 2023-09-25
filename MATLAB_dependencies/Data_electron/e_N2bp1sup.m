function cross_section = e_N2bp1sup(Ep)
% Xs = e_N2bp1sup(E)
% 
% e** + N2 -> e* + N2(b'1Su+) electron excitation cross section (m^2)
% E electron energy (eV)

for iE = numel(Ep):-1:1,
  if Ep(iE) > 12.85 & Ep(iE) < 152,
    cross_section(iE) = (1-12.85/Ep(iE))*exp(-42.893187+3.656131*log(Ep(iE))-0.9243795*log(Ep(iE))^2+0.0642402*log(Ep(iE))^3);
  elseif Ep(iE) >= 152,
    cross_section(iE) = 3.649e-16*log(0.0573*Ep(iE))/Ep(iE);
  else
    cross_section(iE) = 0;
  end
end

cross_section = cross_section*1e-4;