function ThatWentOK = Ie_ztE2Q_zt(nN2,nO2,nO,h_atm, ne, Te)
% IE_ZTE2Q_ZT - Calculate the volume emission/excitation rates
%   Ie_ztE2Q_zt loads all electron-fluxes in a directory and
%   calculates volume emission and excitation-rates for the auroral
%   emissions at 4278 AA, 6730 AA, 7774 AA, 8446 AA, and the O1D, O1S
%   and the N2A3 excitation rates. The function automatically
%   handles the calculations of excitation cross-sections etc.
% 
% Calling:
%  ThatWentOK = Ie_ztE2Q_zt(nN2,nO2,nO,h_atm)
% Input:
%  nN2 - molecular nitrogen number-density (/m^3) profiles [nZ x 1]
%  nO2 - molecular oxygen number-density (/m^3) profiles [nZ x 1]
%  nO  - atomic oxygen number-density (/m^3) profiles [nZ x 1]
%  h_atm - altitude profils [m], [nZ x 1]
% Output
%  ThatWentOK - bolean output, returns 1 if processing worked out

%   Copyright  2018-2019 Bjorn Gustavsson, <bjorn.gustavsson@uit.no>
%   This is free software, licensed under GNU GPL version 2 or later



try
  [t,h_atm,E,mu_lims,Ie_ZTE] = Ie_ztE_loader({'.'});
  szIzte = size(Ie_ZTE);
  % load neutral_atm.mat
  
  emXS4278 = exc_4278(E);
  emXS6730 = exc_6730_N2(E);
  emXS8446_O = exc_8446_O(E);
  emXS8446_O2 = exc_8446_O2(E);
  emX7774_O = exc_7774_O(E);
  emX7774_O2 = exc_7774_O2(E);
  excXS_O1D = exc_O1D(E);
  excXS_O1S = exc_O1S(E);
  dE = diff(E);
  dE = dE([1:end,end]);
  XsO  = get_all_xs('O',E+dE/2);
  XsO2 = get_all_xs('O2',E+dE/2);
  XsN2 = get_all_xs('N2',E+dE/2);
  load N2_levels.dat
  load O2_levels.dat
  load O_levels.dat
  XsOi = O_levels(:,2)'*XsO;
  XsO2i = O2_levels(:,2)'*XsO2;
  XsN2i = N2_levels(:,2)'*XsN2;
  XsN2A3 = XsN2(13,:);

  
  QO = zeros(size(XsO, 1), size(h_atm, 1));
  for ixso = 1:size(XsO, 1)
      Xs = XsO(ixso, :);
      Q = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO,Xs(1:szIzte(3)));
      QO(ixso, :) = Q;
  end

  QO2 = zeros(size(XsO2, 1), size(h_atm, 1));
  for ixso2 = 1:size(XsO2, 1)
      Xs = XsO2(ixso2, :);
      Q = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO2,Xs(1:szIzte(3)));
      QO2(ixso2, :) = Q;
  end

  QN2 = zeros(size(XsN2, 1), size(h_atm, 1));
  for ixsn2 = 1:size(XsN2, 1)
      Xs = XsN2(ixsn2, :);
      Q = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nN2,Xs(1:szIzte(3)));
      QN2(ixsn2, :) = Q;
  end

  save Q_major.mat QO QO2 QN2
  
  Q4278 = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nN2,emXS4278(1:szIzte(3)));
  
  Q6730 = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nN2,emXS6730(1:szIzte(3)));
  
  Q8446 = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO,emXS8446_O(1:szIzte(3))) + ...
          exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO2,emXS8446_O2(1:szIzte(3)));
  Q7774 = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO,emX7774_O(1:szIzte(3))) + ...
          exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO2,emX7774_O2(1:szIzte(3)));
  
  Q8446_O = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO,emXS8446_O(1:szIzte(3)));
  Q7774_O = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO,emX7774_O(1:szIzte(3)));

  Q8446_O2 = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO2,emXS8446_O2(1:szIzte(3)));
  Q7774_O2 = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO2,emX7774_O2(1:szIzte(3)));
  
  % This is to account for quenching for O1D emission %
  % In case Tn was not saved, you must save the msis_file as a .dat and change the name string below
  % msis_file = "/home/etienne/Documents/Julia/AURORA.jl/internal_data/data_neutrals/msis_20051008-2200_70N-19E.dat";
  % load(msis_file);
  % OPS.atmosphere = interp1(msis_20051008_2200_70N_19E(:,1),msis_20051008_2200_70N_19E,h_atm/1e3,'pchip');
  % Tn = OPS.atmosphere(:,end);

  % tauO1D = tau_O1D(nO / 1e6, nO2 / 1e6, nN2 / 1e6, ne / 1e6, Tn, Te, h_atm);
  % QO1D = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO,excXS_O1D(1:szIzte(3))) .* tauO1D;
  % === %
  QO1D = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO,excXS_O1D(1:szIzte(3)));
  QO1S = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO,excXS_O1S(1:szIzte(3)));

  QN2i = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nN2,XsN2i(1:szIzte(3)));
  QO2i = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO2,XsO2i(1:szIzte(3)));
  QOi = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nO,XsOi(1:szIzte(3)));
  
  QN2A3 = exc_tz_of_Ie_ztE(t,h_atm,E(1:szIzte(3)),Ie_ZTE,nN2,XsN2A3(1:szIzte(3)));
  
  save Qzt_all_L.mat t h_atm Q4278 Q6730 Q8446 Q7774 QO1D QO1S QN2i QN2A3 QO2i QOi Q8446_O Q7774_O Q8446_O2 Q7774_O2
  
  c_o_mu = mu_avg(mu_lims);
  idx_z = numel(h_atm);
  [J_down,J_up,IeEdown,IeEup] = Ie2currents(t,h_atm,E(1:szIzte(3)),Ie_ZTE,c_o_mu,idx_z);
  save J.mat t J_down J_up IeEdown IeEup
  
  ThatWentOK = 1;

catch
  ThatWentOK = 0;
end
