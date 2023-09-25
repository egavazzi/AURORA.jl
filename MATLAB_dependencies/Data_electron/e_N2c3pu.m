function Xs = e_N2c3pu(Ep)
% Xs = e_N2c3pu(E)
% 
% e_N2c3pu - electron excitation cross section (m^2)
% E electron energy (eV)
% Cross-sections for energies below 50 eV digitized from Itikawa et
% al 1986, figure 7.6, data from Cartwright et al. (1977),
% renormalized by Trajmar et al. (1983). /BG 20180603

for ie=1:length(Ep),
  if Ep(ie) > 11.03 & Ep(ie) < 20.169,
    cross_section(ie)=(1-11.03/Ep(ie))*exp(-9134.460+12303.41*log(Ep(ie))-6226.840*log(Ep(ie))^2+1397.931*log(Ep(ie))^3-117.4893*log(Ep(ie))^4);
  elseif Ep(ie) >=20.169 & Ep(ie) < 100,
    cross_section(ie)=exp(-145.6415+117.1985*log(Ep(ie))-47.66838*log(Ep(ie))^2+8.547379*log(Ep(ie))^3-0.5779936*log(Ep(ie))^4);    
  elseif Ep(ie) >=100,
    cross_section(ie)=5.5409940e-13/Ep(ie)^3;
  else
    cross_section(ie)=0;
  end
end


E = [12.1 12.67 13.33 13.5 15   17   20   30   50];
s = [ 0.5  1     2     3    3.75 2.05 1.44 0.54 0.21]*1e-21;%-17;

Xs = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
I = find(~isfinite(Xs));
Xs(I) = 0;

inans = find(Ep<11.032);
Xs(inans) = 0;

Xs(Ep>24.4) = cross_section(Ep>24.4)/1e4;
