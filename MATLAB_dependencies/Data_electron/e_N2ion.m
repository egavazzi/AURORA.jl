function cross_section = e_N2ion(Energy)
% Xs = e_N2ion(E)
% 
% e_N2ion - electron ionisation cross section (m^2)
% E electron energy (eV)


for ie=1:length(Energy),
  if Energy(ie) > 15.58 & Energy(ie) < 42.71,
    cross_section(ie)=(1-15.58/Energy(ie))*exp(-833.9598+879.6264*log(Energy(ie))-363.6978*log(Energy(ie))^2+66.81782*log(Energy(ie))^3-4.600032*log(Energy(ie))^4);
  elseif Energy(ie) >=42.71 & Energy(ie) < 250,
    cross_section(ie)=exp(-99.91206+51.41718*log(Energy(ie))-15.61180*log(Energy(ie))^2+2.119856*log(Energy(ie))^3-0.1089437*log(Energy(ie))^4);
  elseif Energy(ie) >=250,
    cross_section(ie)=2.050132e-14*log(0.0275*Energy(ie))/Energy(ie);
  else
    cross_section(ie)=0;
  end
end
cross_section = cross_section/1e4;

% $$$ 
% $$$ E = [17.5 22.5 27.5 35   45   55   65   75  100  150  200  300 400  500  600 700  800   900  1000];
% $$$ s = [2.6   6.4 10.4 14.3 16.8 18.1 18.9 19.3 19.8 18.6 17.4 15  13.1 11.2 10   9.3  8.4   7.9   7.4]*1e-21;%-17;
% $$$ 
% $$$ 
% $$$ Xs = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
% $$$ I = find(~isfinite(Xs));
% $$$ Xs(I) = 0;
% $$$ 
% $$$ Xs = Xs;
% $$$ 
% $$$ inans = find(Ep<15.58);
% $$$ Xs(inans) = 0;

