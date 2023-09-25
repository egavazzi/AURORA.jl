function cs = e_Oelastic(Energy)
% Xs = e_Oelastic(E)
% 
% e_Oelastic - elastic electron collision cross section (m^2)
% E electron energy (eV)

% $$$ xs =[4.1461       6.1685       7.0787       7.2809       7.2809       7.0787]*1e-20;
% $$$ E_o = [1 2 4 6 8 10];
% $$$ 
% $$$ Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
% $$$ i = find(~isfinite(Xs));
% $$$ Xs(i) = 0;
% $$$ i = find(E<.01);
% $$$ Xs(i) = 0;

for ie=1:length(Energy),
  if Energy(ie) > 0.2 & Energy(ie) < 10000,
    cs(ie)=2.8e-17*exp(2.8+0.39*log(Energy(ie))-0.078*log(Energy(ie))^2-0.0102*log(Energy(ie))^3+0.00095*log(Energy(ie))^4); 
  elseif Energy(ie) >= 10000.,
    cs(ie)=2.8e-17*2572.3/Energy(ie);
  end
end    
cs = cs/1e4;
