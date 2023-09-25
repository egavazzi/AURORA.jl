function cross_section = e_N2elastic(Energy)
% cross_section =  e_N2elastic(Energy)
% 
% e_N2elastic - elastic electron collision cross section (m^2)
% E electron energy (eV)


for ie=1:length(Energy),
  if Energy(ie) >= 1 & Energy(ie) < 1.733,
    cross_section(ie)=exp(-34.56132+0.9475102*log(Energy(ie))-4.235151*log(Energy(ie))^2+6.523666*log(Energy(ie))^3);
  elseif Energy(ie) >= 1.733 & Energy(ie) < 2.006,
    cross_section(ie)=exp(9.113062-242.3374*log(Energy(ie))+441.5981*log(Energy(ie))^2-262.3439*log(Energy(ie))^3);
  elseif Energy(ie) >= 2.006 & Energy(ie) < 2.317,
    cross_section(ie)=exp(72.26128-455.1942*log(Energy(ie))+642.8157*log(Energy(ie))^2-299.3418*log(Energy(ie))^3);
  elseif Energy(ie) >= 2.317 & Energy(ie) < 2.554,
    cross_section(ie)=exp(1767.169-6203.92*log(Energy(ie))+7114.73*log(Energy(ie))^2-2716.333*log(Energy(ie))^3);
  elseif Energy(ie) >= 2.554 & Energy(ie) < 2.775,
    cross_section(ie)=exp(867.3897-2850.537*log(Energy(ie))+3000.209*log(Energy(ie))^2-1050.929*log(Energy(ie))^3);
  elseif Energy(ie) >= 2.775 & Energy(ie) < 2.981,
    cross_section(ie)=exp(679.3589-2066.956*log(Energy(ie))+1995.594*log(Energy(ie))^2-641.9862*log(Energy(ie))^3);
  elseif Energy(ie) >= 2.981 & Energy(ie) < 3.215,
    cross_section(ie)=exp(-400.2765+865.4071*log(Energy(ie))-668.1613*log(Energy(ie))^2+167.3741*log(Energy(ie))^3);
  elseif Energy(ie) >= 3.215 & Energy(ie) < 3.457,
    cross_section(ie)=exp(1815.893-4716.509*log(Energy(ie))+4004.182*log(Energy(ie))^2-1132.109*log(Energy(ie))^3);
  elseif Energy(ie) >= 3.457 & Energy(ie) < 4.33,
    cross_section(ie)=exp(45.31069-184.5681*log(Energy(ie))+142.8378*log(Energy(ie))^2-36.88586*log(Energy(ie))^3);
  elseif Energy(ie) >= 4.33 & Energy(ie) < 1000,
    cross_section(ie)=exp(-35.26294+0.8019902*log(Energy(ie))-0.202546*log(Energy(ie))^2+0.007484*log(Energy(ie))^3);
  elseif Energy(ie) >= 1000,
    cross_section(ie)=9.234e-14/Energy(ie);
  else
    cross_section(ie)=0;
  end
end
cross_section = cross_section/1e4;
% $$$ 
% $$$ xs =[  2.1444   4.2292   5.5908   8.1113   10.235]*1e-20;
% $$$ E_o = [0.01     0.05     0.1      0.3       0.9];
% $$$ xs = [xs [10.043 11.466 19.224 16.315 18.578 22.522 18.319 25.172 16.832 20.517 17.026 17.672 15.927 16.897 12.823 15.022 12.629 13.405 12.629 12.888]*1e-20];
% $$$ E_o =[E_o [1      1.5617 1.9191 1.9957 2.0723 2.2    2.3277 2.4119 2.583  2.6723 2.8    2.8766 2.9532 3.0681 3.2213 3.3489 3.4638 3.6043 3.783 3.8851]];
% $$$ xs = [xs [11.444   12.101   9.7701   6.7342   5.5908   2.6561   1.0476]*1e-20];
% $$$ E_o =[E_o [5       10      30       70      100      300     1000]];
% $$$ Xs = exp([interp1(log(E_o),log(xs),log(E(E<E_o(end))),'pchip') interp1(log(E_o),log(xs),log(E(E>=E_o(end))),'linear','extrap')]);
% $$$ i = find(~isfinite(Xs));
% $$$ Xs(i) = 0;
% $$$ i = find(E<1e-2);
% $$$ Xs(i) = 0;
