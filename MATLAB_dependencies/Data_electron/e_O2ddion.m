function cross_section = e_O2ddion(Ep)
% Xs = e_O2dion(E)
% 
% e_O2dion - dissociative-double ionization cross section (m^2)
% E electron energy (eV)


for iE = numel(Ep):-1:1
  if Ep(iE) > 32.51 & Ep(iE) < 300,
    cross_section(iE) = (1-32.51/Ep(iE))*exp(-228.0913+141.5861*log(Ep(iE))-40.00023*log(Ep(iE))^2+5.079899*log(Ep(iE))^3-0.2454952*log(Ep(iE))^4);
  elseif Ep(iE) >= 200,
    cross_section(iE) = 8.141987e-16*log(0.12511*Ep(iE))/Ep(iE);
  else
    cross_section(iE) = 0;
  end
end
cross_section = cross_section/1e4;