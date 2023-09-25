function cross_section = e_N2e3sgp(Ep)
% Xs = e_N2e3sgp(E)
% 
% e_N2e3sgp - electron excitation cross section (m^2)
% E electron energy (eV)


for ie=1:length(Ep),
  if Ep(ie) > 11.88 & Ep(ie) < 30.35,
    cross_section(ie)=(1-11.88/Ep(ie))*exp(74.2133-124.664*log(Ep(ie))+44.6708*log(Ep(ie))^2-5.30151*log(Ep(ie))^3);
  elseif Ep(ie) >=30.35 & Ep(ie) < 50,
    cross_section(ie)=exp(-105.2329+54.4886*log(Ep(ie))-14.77986*log(Ep(ie))^2+1.23869*log(Ep(ie))^3);    
  elseif Ep(ie) >=50,
    cross_section(ie)=8.7619e-15/Ep(ie)^3;
  else
    cross_section(ie)=0;
  end
end
cross_section = cross_section/1e4;

% $$$ E = [12.255 15  17  20  30 50];
% $$$ s = [.01    .2  .2  .5  .4 .1]*1e-22;%-18;
% $$$ 
% $$$ Xs = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
% $$$ I = find(~isfinite(Xs));
% $$$ Xs(I) = 0;
% $$$ 
% $$$ inans = find(Ep<11.875);
% $$$ Xs(inans) = 0;
% $$$ 
