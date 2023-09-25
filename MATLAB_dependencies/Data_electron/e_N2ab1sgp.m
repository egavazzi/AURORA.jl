function cross_section = e_N2ab1sgp(Ep)
% Xs =  e_N2ab1sgp(E)
% 
% e_N2ab1sgp - electron excitation cross section (m^2)
% E electron energy (eV)


for ie=1:length(Ep),
  if Ep(ie) > 12.25 & Ep(ie) < 24.98,
    cross_section(ie)=(1-12.25/Ep(ie))*exp(91.77961-148.1616*log(Ep(ie))+55.75255*log(Ep(ie))^2-6.95604*log(Ep(ie))^3);
  elseif Ep(ie) >=24.98 & Ep(ie) < 50,
    cross_section(ie)=exp(80.2784-92.0627*log(Ep(ie))+23.36969*log(Ep(ie))^2-1.985404*log(Ep(ie))^3);    
  elseif Ep(ie) >=50,
    cross_section(ie)=7.143e-17/Ep(ie);
  else
    cross_section(ie)=0;
  end
end
cross_section = cross_section/1e4;

I = find(~isfinite(cross_section));
cross_section(I) = 0;

% $$$ E = [12.255 15  17  20  30 50];
% $$$ s = [.01    3.8 3.7 5.8  2  1.3]*1e-22;%-18;
% $$$ 
% $$$ Xs = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
% $$$ I = find(~isfinite(Xs));
% $$$ Xs(I) = 0;
% $$$ 
% $$$ inans = find(Ep<12.255);
% $$$ Xs(inans) = 0;
