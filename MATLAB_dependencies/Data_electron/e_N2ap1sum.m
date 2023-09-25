function Xs = e_N2ap1sum(Ep)
% Xs = e_N2ap1sum(E)
% 
% e_N2ap1sum - electron excitation cross section (m^2)
% E electron energy (eV)
% Cross-sections for energies below 50 eV digitized from Itikawa et
% al 1986, figure 7.12, data from Cartwright et al. (1977),
% renormalized by Trajmar et al. (1983). /BG 20180603

for ie=1:length(Ep),
  if Ep(ie) > 8.4 & Ep(ie) < 29.937,
    cross_section(ie)=(1-8.4/Ep(ie))*exp(-2652.316+3473.205*log(Ep(ie))-1723.654*log(Ep(ie))^2+378.6041*log(Ep(ie))^3-31.07141*log(Ep(ie))^4);
  elseif Ep(ie) >=29.937 & Ep(ie) < 100,
    cross_section(ie)=exp(-436.8312+412.7602*log(Ep(ie))-159.5703*log(Ep(ie))^2+27.18180*log(Ep(ie))^3-1.728542*log(Ep(ie))^4);    
  elseif Ep(ie) >=100,
    cross_section(ie)=2.495837e-17/Ep(ie);
  else
    cross_section(ie)=0;
  end
end

E = [8.3987 15  18  20    30 50];
s = [.0001  11   5   3.5   2  .8]*1e-22;%-18;

Xs = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
I = find(~isfinite(Xs));
Xs(I) = 0;

inans = find(Ep<8.3987);
Xs(inans) = 0;

Xs(Ep>45) = cross_section(Ep>45)/1e4;
