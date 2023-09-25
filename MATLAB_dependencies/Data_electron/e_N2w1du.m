function Xs = e_N2w1du(Ep)
% Xs = e_N2w1du(E)
% 
% e_N2w1du - electron excitation cross section (m^2)
% E electron energy (eV)
% Cross-sections for energies below 50 eV digitized from Itikawa et
% al 1986, figure 7.12, data from Cartwright et al. (1977),
% renormalized by Trajmar et al. (1983). /BG 20180603

for ie=1:length(Ep),
  if Ep(ie) > 8.89 & Ep(ie) < 20.487,
    cross_section(ie)=(1-8.89/Ep(ie))*exp(-5231.492+7031.762*log(Ep(ie))-3566.315*log(Ep(ie))^2+802.9877*log(Ep(ie))^3-67.74330*log(Ep(ie))^4);
  elseif Ep(ie) >=20.487 & Ep(ie) < 100,
    cross_section(ie)=exp(-131.5858+104.4165*log(Ep(ie))-43.27655*log(Ep(ie))^2+7.735902*log(Ep(ie))^3-0.5085983*log(Ep(ie))^4);    
  elseif Ep(ie) >=100,
    cross_section(ie)=7.415168e-17/Ep(ie);
  else
    cross_section(ie)=0;
  end
end


E = [8.8895 12  15  17  20  30    50];
s = [.0001  12  9.5  7   3   1.75  0.7]*1e-22;%-18;

Xs = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
I = find(~isfinite(Xs));
Xs(I) = 0;

inans = find(Ep<8.8895);
Xs(inans) = 0;

Xs(30<=Ep&Ep<35) = ( Xs(30<=Ep&Ep<35)*2/3 + ...
                    cross_section(30<=Ep&Ep<35)*1/3/1e4 );
Xs(35<=Ep&Ep<40) = ( Xs(35<=Ep&Ep<40)*2/3 + ...
                    cross_section(35<=Ep&Ep<40)*1/3/1e4 );
Xs(40<=Ep) = cross_section(40<=Ep)/1e4;
