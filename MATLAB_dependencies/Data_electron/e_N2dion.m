function cross_section = e_N2dion(Ep)
% Xs = e_N2dion(E)
% 
% e_N2dion - dissociative ionization cross section (m^2)
% E electron energy (eV)

% $$$ E = [30  40   50  60  70  80  90 100 110 120 140 160 180 200 250 300];
% $$$ s = [0.16 0.97 2.5 3.9 4.6 5.0 5.2 5.4 5.5 5.5 5.2 4.8 4.6 4.4 3.9 3.4]*1e-21;%-17;
% $$$ 
% $$$ 
% $$$ Xs = exp([interp1(E,log(s),Ep(Ep<E(end)),'pchip') interp1(log(E),log(s),log(Ep(Ep>=E(end))),'linear','extrap')] );
% $$$ I = find(~isfinite(Xs));
% $$$ Xs(I) = 0;
% $$$ 
% $$$ Xs = Xs;
% $$$ 
% $$$ inans = find(Ep<25);
% $$$ Xs(inans) = 0;

for iE = numel(Ep):-1:1,
  if Ep(iE) > 24.00 & Ep(iE) < 42.82,
    cross_section(iE) = (1-24.00/Ep(iE))*exp(2066.673-2216.257*log(Ep(iE))+870.8756*log(Ep(iE))^2-151.3981*log(Ep(iE))^3+9.830359*log(Ep(iE))^4);
  elseif Ep(iE) >= 42.82 & Ep(iE) < 600,
    cross_section(iE) = exp(-159.5311+90.55610*log(Ep(iE))-25.06691*log(Ep(iE))^2+3.076254*log(Ep(iE))^3-0.1420719*log(Ep(iE))^4);
  elseif Ep(iE) >= 600,
    cross_section(iE) = 3.526978e-15*log(0.0275*Ep(iE))/Ep(iE);
  else
    cross_section(iE) = 0;
  end
end

cross_section = cross_section*1e-4;