function cross_section = e_O2ionb4sgm(Energy)
% Xs = e_O2ionb4sgm(E)
% 
% e_O2ionb4sgm - electron ionisation cross section (m^2) to the
% excited states O_2^+(b^4\Sigma_g^-)
% E electron energy (eV)


for iE = numel(Energy):-1:1
  if Energy(iE) > 18.2 & Energy(iE) < 250,
    cross_section(iE)=(1-18.2/Energy(iE))*exp(-26.45362-12.67161*log(Energy(iE))+4.799453*log(Energy(iE))^2-0.7680145*log(Energy(iE))^3+0.04377051*log(Energy(iE))^4);
  elseif Energy(iE) >= 250,
    cross_section(iE)=1.980794e-15*log(0.024825*Energy(iE))/Energy(iE);
  else
    cross_section(iE)=0;
  end
end
    
cross_section = cross_section/1e4;
