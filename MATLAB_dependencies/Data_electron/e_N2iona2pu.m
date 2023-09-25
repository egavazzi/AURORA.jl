function cross_section = e_N2iona2pu(Energy)
% Xs = e_N2iona2pu(E)
% 
% e_N2ionA2Pu - electron ionisation cross section (m^2) to the
% first electronically excited state of N2+
% E electron energy (eV)


for iE = numel(Energy):-1:1,
  if Energy(iE) > 16.73 & Energy(iE) < 42.85,
    cross_section(iE) = (1-16.73/Energy(iE))*exp( -40.81558-6.435371*log(Energy(iE))+6.032905*log(Energy(iE))^2-1.545984*log(Energy(iE))^3+0.1261087*log(Energy(iE))^4);
  elseif Energy(iE) >= 42.85 & Energy(iE) < 300,
    cross_section(iE) = exp(-63.63995+19.39587*log(Energy(iE))-5.333680*log(Energy(iE))^2+0.6665464*log(Energy(iE))^3-0.03254920*log(Energy(iE))^4);
  elseif Energy(iE) >= 300,
    cross_section(iE) = 9.022592e-15*log(0.0275*Energy(iE))/Energy(iE);
  else
    cross_section(iE) = 0;
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

