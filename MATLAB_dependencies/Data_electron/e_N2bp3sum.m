function Xs = e_N2bp3sum(Ep)
% Xs = e_N2bp3sum(E)
% 
% e_N2bp3sum - electron excitation cross section (m^2)
% E electron energy (eV)
% Cross-sections for energies below 50 eV digitized from Itikawa et
% al 1986, figure 7.12, data from Cartwright et al. (1977),
% renormalized by Trajmar et al. (1983). /BG 20180603


for ie=1:length(Ep),
  if Ep(ie) > 8.16 & Ep(ie) < 29.55,
    cross_section(ie)=(1-8.16/Ep(ie))*exp(-667.3893+764.3953*log(Ep(ie))-346.0872*log(Ep(ie))^2+69.19737*log(Ep(ie))^3-5.167390*log(Ep(ie))^4);
  elseif Ep(ie) >=29.55 & Ep(ie) < 100,
    cross_section(ie)=exp(-34.27635+0.0*log(Ep(ie))-1.745869*log(Ep(ie))^2+0.5473512*log(Ep(ie))^3-0.0553294*log(Ep(ie))^4);    
  elseif Ep(ie) >=100,
    cross_section(ie)=2.77e-13/Ep(ie)^3;
  else
    cross_section(ie)=0;
  end
end

E = [8.1647 12 15  18  20  30 50];
s = [.0001  8  13   8   4   3  2]*1e-22;%-18;

Xs = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
I = find(~isfinite(Xs));
Xs(I) = 0;

inans = find(Ep<8.1647);
Xs(inans) = 0;

Xs(Ep>30) = cross_section(Ep>30)/1e4;
