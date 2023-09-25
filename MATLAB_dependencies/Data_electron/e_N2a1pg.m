function Xs = e_N2a1pg(Ep)
% Xs = e_N2a1pg(E)
% 
% e_N2a1pg - electron excitation cross section (m^2)
% E electron energy (eV)
% Cross-sections for energies below 50 eV digitized from Itikawa et
% al. (1986), figure 7.4, data from Finn and Doering (1976), and
% Cartwright et al. (1977), renormalized by Trajmar et
% al. (1983). /BG 20180603 

%  Copyright © Bjorn Gustavsson 20100527, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

for ie=1:length(Ep),
  if Ep(ie) > 8.55 & Ep(ie) < 100.0,
    cross_section(ie)=(1-8.55/Ep(ie))*exp(-108.6546+73.07788*log(Ep(ie))-27.46333*log(Ep(ie))^2+4.465812*log(Ep(ie))^3-0.2689957*log(Ep(ie))^4);
  elseif Ep(ie) >=100.0,
    cross_section(ie)=7.207318e-16/Ep(ie);
  else
    cross_section(ie)=0;
  end
end

E = [10  12.5 15   16.25 17.5 20   22.0 30.0 35.0 40.0 50    60   75   90  100];
s = [0.42 2.29 3.67 4     3.67 3.08 2.67 2.12 1.67 1.5  1.17  0.83 0.75 0.57 0.5]*1e-21;%-17

Xs = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
I = find(~isfinite(Xs));
Xs(I) = 0;

inans = find(Ep<8.5489);
Xs(inans) = 0;

Xs(Ep>31.25) = cross_section(Ep>31.25)/1e4;
