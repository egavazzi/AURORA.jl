function cross_section = e_N2ddion(Ep)
% Xs = e_N2ddion(E)
% 
% e_N2ddion - double dissociative ionization cross section (m^2)
% i.e. N2 + e** -> e* + N^+ + N^+
% E electron energy (eV)


for iE = numel(Ep):-1:1,
  if Ep(iE) > 42.00 & Ep(iE) < 44.45,
    cross_section(iE) = (1-42.00/Ep(iE))*exp(830.9748-744.6075*log(Ep(iE))+235.7621*log(Ep(iE))^2-32.65055*log(Ep(iE))^3+1.664205*log(Ep(iE))^4);
  elseif Ep(iE) >= 44.45 & Ep(iE) < 550,
    cross_section(iE) = exp(-240.8384+137.2195*log(Ep(iE))-34.66862*log(Ep(iE))^2+3.873074*log(Ep(iE))^3-0.1624777*log(Ep(iE))^4);
  elseif Ep(iE) >= 550,
    cross_section(iE) = 9.873940e-16*log(0.0275*Ep(iE))/Ep(iE);
  else
    cross_section(iE) = 0;
  end
end    

cross_section = cross_section*1e-4;