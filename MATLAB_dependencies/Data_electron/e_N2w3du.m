function Xs = e_N2w3du(Ep)
% Xs = e_N2w3du(E)
% 
% e_N2w3du - electron excitation cross section (m^2)
% E electron energy (eV)
% Cross-sections for energies below 50 eV digitized from Itikawa et
% al 1986, figure 7.1, data from Cartwright et al. (1977),
% renormalized by Trajmar et al. (1983). /BG 20180603

%  Copyright © Bjorn Gustavsson 20180603, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later
 
for ie=1:length(Ep),
  if Ep(ie) > 7.36 & Ep(ie) < 100,
    cross_section(ie)=(1-7.36/Ep(ie))*exp(-275.1462+274.0435*log(Ep(ie))-116.3761*log(Ep(ie))^2+21.57467*log(Ep(ie))^3-1.484377*log(Ep(ie))^4);
  elseif Ep(ie) >=100,
    cross_section(ie)=7.036322e-13/Ep(ie)^3;
  else
    cross_section(ie)=0;
  end
end

E = [7  8 9  10   12.5 15 17 20 30 50];
s = [.2 2 4.5 7.5 24.5 33 34 31  7  2.2]*1e-22;%-18;

Xs = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
I = find(~isfinite(Xs));
Xs(I) = 0;

inans = find(Ep<7.3622);
Xs(inans) = 0;

Xs(Ep>28.25) = cross_section(Ep>28.25)/1e4;


  
