function Xs = e_N2a3sup(Ep)
% Xs = e_N2a3sup(E)
% 
% e_N2a3sup - electron excitation cross section (m^2)
% E electron energy (eV)

E = [7  8 9   10 13   15   17 20 30 50];
s = [.2 2 4.5 15 17.6 20.5 22 15  6  3.8]*1e-22;%-18;



for ie=1:length(Ep),
  if Ep(ie) > 6.17 & Ep(ie) < 11.096,
    Xs(ie)=(1-6.17/Ep(ie))*exp(-48.28302-501.0230*log(Ep(ie))+635.6719*log(Ep(ie))^2-264.7979*log(Ep(ie))^3+36.53586*log(Ep(ie))^4);
  elseif Ep(ie) >= 11.096 & Ep(ie) < 17.05,
    Xs(ie)=exp(4550.309-6733.425*log(Ep(ie))+3697.708*log(Ep(ie))^2-900.7303*log(Ep(ie))^3+82.11416*log(Ep(ie))^4);
  elseif Ep(ie) >= 17.05 & Ep(ie) < 100.0,
    Xs(ie)=exp(-142.1066+124.6773*log(Ep(ie))-55.47813*log(Ep(ie))^2+10.80510*log(Ep(ie))^3-0.7832104*log(Ep(ie))^4);
  elseif Ep(ie) >=100,
    Xs(ie)=9.608389e-13/Ep(ie)^3;
  else
    Xs(ie)=0;
  end
end

XsItikawa = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
I = find(~isfinite(XsItikawa));
XsItikawa(I) = 0;

Xs = [XsItikawa(Ep<=25),Xs(Ep>25)/1e4];
inans = find(Ep<6.1688);
Xs(inans) = 0;
