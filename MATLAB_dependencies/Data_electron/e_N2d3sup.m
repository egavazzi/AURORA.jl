function cross_section = e_N2d3sup(Ep)
% Xs = e_N2d3sup(E)
% 
% e_N2D3Su+ - electron excitation cross section (m^2)
% E electron energy (eV)

for iE = numel(Ep):-1:1,
  if Ep(iE) > 12.85 & Ep(iE) < 50,
    cross_section(iE) = (1-12.85/Ep(iE))*exp(0.938347-31.08899*log(Ep(iE))+8.257418*log(Ep(iE))^2-0.793736*log(Ep(iE))^3);
  elseif Ep(iE) >= 50,
    cross_section(iE) = 6.3166e-14/Ep(iE)^3;
  else
    cross_section(iE) = 0;
  end
end

cross_section = cross_section*1e-4;
