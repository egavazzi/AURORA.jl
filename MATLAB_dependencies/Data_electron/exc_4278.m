function sigma = exc_4278(E)
% Xs = exc_4278(E)
% 
% exc_4278 - emission cross section (m^2) for 4278 AA
% from molecular nitrogen ions. E electron energy (eV)

% $$$ xs_n2p1n = [.1 2.1 5.5 3]*1e-22;%-18;
% $$$ E_n2p1n = [18 30 90 400];
% $$$ X_n2p1n = exp(interp1(log(E_n2p1n),log(xs_n2p1n),log(E),'cubic'));
% $$$ i = find(~isfinite(X_n2p1n));
% $$$ X_n2p1n(i) = 0;
% $$$ inans = find(E<18.75);
% $$$ X_n2p1n(inans) = 0;
% $$$ 
% $$$ sigma = X_n2p1n;


for ie=1:length(E),
  if E(ie) > 18.75 & E(ie) < 43.4,
    cross_section(ie)=(1-18.75/E(ie))*exp(238.9714-290.0367*log(E(ie))+112.9749*log(E(ie))^2-19.41400*log(E(ie))^3+1.241946*log(E(ie))^4);
  elseif E(ie) >=43.4 & E(ie) < 300,
    cross_section(ie)=exp(-65.15851+19.39587*log(E(ie))-5.333680*log(E(ie))^2+0.6665464*log(E(ie))^3-0.03254920*log(E(ie))^4);
  elseif E(ie) >=300,
    cross_section(ie)=1.976205e-15*log(0.0275*E(ie))/E(ie);
  else
    cross_section(ie)=0;
  end
end

sigma = cross_section/1e4/3.3773;
